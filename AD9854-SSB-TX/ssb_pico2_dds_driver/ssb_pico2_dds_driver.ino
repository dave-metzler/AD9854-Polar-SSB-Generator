// SSB polar modulator -- AD9854 DRIVER STAGE v3, for Raspberry Pi Pico 2 (RP2350).
// UNTESTED on hardware (design verified in a Python model; see notes at the bottom of this header).
// Part 2 of 2: receives baseband I/Q samples (44.1 kHz) from the Teensy 4.1 DSP stage
// (ssb_teensy41_audio_dsp.ino) over UART, converts them to polar form HERE, and drives the AD9854
// with a fixed carrier frequency plus a per-update PHASE-OFFSET word and AMPLITUDE word.
//
// WHY THIS REWRITE (v3) -- what was wrong with v2
// ------------------------------------------------
// Simulating the v2 chain with speech-like input gave only ~25 dB opposite-sideband rejection and a
// carrier ~30 dB below the wanted sideband, while the same Hilbert front end alone is good for ~60 dB.
// The loss was in the polar/DDS back end:
//  1. UPDATE RATE.  The DDS was updated at 8.82 kHz.  A polar modulator is a sampled system: holding
//     (envelope, phase) for 113 us puts sampling images at +-(8820 - f_audio) around the carrier --
//     for a 3 kHz component that is only ~15 dB down and lands 5.8 kHz from the carrier, i.e. right
//     inside what an SDR displays as "opposite sideband / garble".  Envelope/phase are also much wider
//     in bandwidth than the audio they come from, so decimating them to 8.82 kHz aliased real content
//     into the passband.  Now: one update per audio sample (44.1 kHz) x INTERP (default 2 -> 88.2 kHz).
//  2. FREQUENCY-INTEGRATION CONTROL.  The DDS phase is the running integral of the frequency words, so
//     every approximation upstream (gate/hold at envelope nulls, one-pole smoothing, envelope-weighted
//     decimation average, 100 Hz/cycle slew limiter, +-16 kHz clamp) permanently adds phase error that
//     never washes out, and an SSB envelope null needs a ~pi phase jump that frequency steps cannot
//     deliver.  Now the DDS FTW is fixed at the carrier and the AD9854's 14-bit PHASE ADJUST register
//     (0x00/0x01) is written with the absolute phase every update.  No integration, no error build-up,
//     nulls are exact, so the gate, smoother, weighted average and slew limiter are all gone.
//  3. UPDATE INSTANT JITTER.  Register writes are now issued FIRST in each tick from values computed
//     during the previous tick, so the UD edge lands at a fixed delay after the timer interrupt
//     instead of after a data-dependent amount of math.
//  4. RATE LOCK.  The Teensy audio clock and this chip's clock are not the same oscillator.  Updates
//     are paced by a PWM-slice wrap interrupt whose period is dithered by one clk_sys count (6.7 ns)
//     and steered by the (low-passed) fill level of the receive FIFO, so the update rate follows the
//     Teensy audio clock exactly -- no dropped or repeated samples.  (v2 also had a pacing bug on the
//     Teensy side that slowly overflowed its ring buffer.)
//
// AD9854 wiring: UNCHANGED --
//   GP0-7 D0-D7, GP8-13 A0-A5, GP14 WR, GP15 RD, GP16 RST, GP17 S/P,
//   GP18 -> I/O UD CLK (through ~1k series resistor, as before).
//   UART1: GP21 = RX <- Teensy Serial1 TX (pin 1).  GP20 = TX (unused, one-way link).
//   Local PTT: GP22 to GND via the same switch as the Teensy's PTT (redundant hardware mute).
//   Common GND with the Teensy.
//
// LINK FORMAT (4 bytes per sample, self-synchronising, 4 Mbaud): 7 payload bits per byte, byte 0 has
// its MSB set, bytes 1-3 have it clear.  28-bit payload, MSB first:
//     [27] PTT   [26] even-parity bit over payload   [25:13] I (13-bit signed)   [12:0] Q (13-bit signed)
// |I+jQ| = 4095 means full scale.  A damaged byte costs one sample; the receiver resyncs on the next
// byte that has its MSB set.  Raw wire load: 44.1k x 40 bits = 1.76 Mbit/s (44% of 4 Mbaud).
//
// Serial commands (115200 baud, USB):
//   <n>              carrier frequency in MHz, e.g. "14.15"   (default 14.15 MHz)
//   x<Hz>[,<amp>]    bring-up test: replace the Teensy data with a single tone at carrier+<Hz>
//                    (amp 0-4095, default 2000) through the SAME phase/amp/UD path. "x" alone = off.
//                    Expect ONE clean line at carrier+Hz with the carrier and opposite sideband
//                    >50 dB down; if not, the DDS/UD path (not the audio chain) is at fault.
//   i                status
//
// If the ISR cannot keep up (status shows isr max >= period), set INTERP to 1.

#include <Arduino.h>
#include <math.h>
#include "pico/stdlib.h"
#include "pico/multicore.h"
#include "hardware/gpio.h"
#include "hardware/uart.h"
#include "hardware/pwm.h"
#include "hardware/irq.h"
#include "hardware/clocks.h"
#include "hardware/sync.h"

#ifndef PWM_DEFAULT_IRQ_NUM
#define PWM_DEFAULT_IRQ_NUM() PWM_IRQ_WRAP_0
#endif

// ---- Tunables ----
#define INTERP 2                              // DDS updates per input sample (1, 2 or 4). 2 -> 88.2 kHz
const double   SAMPLE_RATE_HZ = 44100.0;      // Teensy 4.x Audio Library rate (44117.647 is the Teensy 3.x value)
const uint32_t LINK_BAUD      = 4000000;
const uint32_t STALE_TIMEOUT_US = 50000;      // no good packet for 50 ms -> mute
const int      RING_TARGET    = 160;          // desired FIFO fill (samples) ~3.6 ms; Teensy sends in bursts
                                               // of up to 128, so the target must exceed ~70 + margin
const float    PLL_GAIN_Q16   = 6000.0f;      // clk_sys counts (Q16) of period change per sample of fill error

// ---- AD9854 parallel wiring (unchanged) ----
const uint32_t BUS_MASK = 0x3FFF;
const uint8_t PIN_WR = 14, PIN_RD = 15, PIN_RST = 16, PIN_SP = 17;
const uint8_t PIN_UD = 18;
const uint8_t PIN_UART_TX = 20, PIN_UART_RX = 21;
const uint8_t PIN_LOCAL_PTT = 22;
const uint PWM_SLICE = 11;                    // free-running slice used only as a timebase

// ---- Clocking ----
const double F_REF = 30e6;
const int    MULT  = 10;
const double F_SYS = F_REF * MULT;            // 300 MHz

// ---- Receive FIFO: core1 (UART) produces, core0 (ISR) consumes ----
#define RING_SIZE 512
#define RING_MASK (RING_SIZE - 1)
static uint32_t s_ring[RING_SIZE];            // raw 28-bit payload words
static volatile uint16_t s_head = 0, s_tail = 0;

volatile uint32_t g_packetCount = 0, g_badPackets = 0, g_ringFull = 0, g_underruns = 0;
volatile uint32_t g_lastPacketUs = 0;

// ---- Shared control (loop() -> ISR) ----
volatile double   g_carrierHz = 14.15e6;
volatile uint32_t g_ftwHi = 0, g_ftwLo = 0;   // pending FTW (48-bit) split in two
volatile bool     g_ftwPending = false;
volatile bool     g_testActive = false;
volatile float    g_testHz = 1000.0f;
volatile uint16_t g_testAmp = 2000;

// ---- Diagnostics ----
volatile float    g_fillAvg = 0;
volatile uint32_t g_isrMaxUs = 0;
volatile uint32_t g_periodCounts = 0;
volatile int16_t  g_lastI = 0, g_lastQ = 0;
volatile bool     g_lastPtt = false;

static uint32_t s_periodNomCounts = 1700;     // clk_sys counts per DDS update, set in setup()
static uint32_t s_periodNomUs = 11;

// ---- AD9854 low-level (same timing as the validated sketches) ----
static inline void __not_in_flash_func(writeReg)(uint8_t addr, uint8_t val) {
  gpio_put_masked(BUS_MASK, ((uint32_t)(addr & 0x3F) << 8) | val);
  busy_wait_at_least_cycles(15);
  gpio_put(PIN_WR, 0);
  busy_wait_at_least_cycles(20);
  gpio_put(PIN_WR, 1);
  busy_wait_at_least_cycles(15);
}

static inline void __not_in_flash_func(writeFTWWords)(uint32_t hi, uint32_t lo) {
  // 48-bit FTW = hi(16) : lo(32), written MSB first to 0x04..0x09
  writeReg(0x04, (hi >> 8) & 0xFF);
  writeReg(0x05, hi & 0xFF);
  writeReg(0x06, (lo >> 24) & 0xFF);
  writeReg(0x07, (lo >> 16) & 0xFF);
  writeReg(0x08, (lo >> 8) & 0xFF);
  writeReg(0x09, lo & 0xFF);
}

static void computeFTW(double fout, uint32_t &hi, uint32_t &lo) {
  uint64_t ftw = (uint64_t)llround(fout / F_SYS * 281474976710656.0);   // 2^48
  hi = (uint32_t)((ftw >> 32) & 0xFFFF);
  lo = (uint32_t)(ftw & 0xFFFFFFFFu);
}

// Phase Adjust Register 1 (0x00 = bits 13:8, 0x01 = bits 7:0) and the two OSK multipliers.
static inline void __not_in_flash_func(writePhaseAmp)(uint16_t ph14, uint16_t amp) {
  writeReg(0x00, (ph14 >> 8) & 0x3F);
  writeReg(0x01, ph14 & 0xFF);
  writeReg(0x21, (amp >> 8) & 0x0F);
  writeReg(0x22, amp & 0xFF);
  writeReg(0x23, (amp >> 8) & 0x0F);
  writeReg(0x24, amp & 0xFF);
}

static inline void __not_in_flash_func(pulseUpdate)() {
  gpio_put(PIN_UD, 1);
  busy_wait_at_least_cycles(20);
  gpio_put(PIN_UD, 0);
  busy_wait_at_least_cycles(15);
}

// atan2 for the ISR: RAM-resident polynomial, max error ~3e-5 rad (14-bit phase LSB is 3.8e-4 rad).
static inline float __not_in_flash_func(fastAtan2)(float y, float x) {
  float ax = fabsf(x), ay = fabsf(y);
  float mx = ax > ay ? ax : ay, mn = ax > ay ? ay : ax;
  if (mx <= 0.0f) return 0.0f;
  float a = mn / mx, s = a * a, c = s * a, q = s * s;
  float r = 0.024840285f * q + 0.18681418f;
  float t = -0.094097948f * q - 0.33213072f;
  r = r * s + t;
  r = r * c + a;
  if (ay > ax) r = 1.57079637f - r;
  if (x < 0.0f) r = 3.14159274f - r;
  if (y < 0.0f) r = -r;
  return r;
}

// ---- ISR state ----
static uint16_t s_nextPh = 0, s_nextAmp = 0;        // register values computed LAST tick, written THIS tick
static float s_z0i = 0, s_z0q = 0, s_z1i = 0, s_z1q = 0;
static int   s_sub = 0;
static bool  s_running = false;                     // FIFO prefilled?
static bool  s_pttBit = false;
static float s_fillAvg = 0.0f;
static uint32_t s_acc = 0;                          // dither accumulator, Q16 counts
static int32_t  s_periodQ16 = 0;
static float s_testPh = 0.0f;
static uint16_t s_lastPh = 0;

static inline void decodeWord(uint32_t w, float &i, float &q, bool &ptt) {
  int32_t iv = (int32_t)((w >> 13) & 0x1FFF);
  int32_t qv = (int32_t)(w & 0x1FFF);
  if (iv & 0x1000) iv -= 0x2000;
  if (qv & 0x1000) qv -= 0x2000;
  i = (float)iv; q = (float)qv;
  ptt = ((w >> 27) & 1) != 0;
}

static void __not_in_flash_func(ddsTick)() {
  pwm_clear_irq(PWM_SLICE);
  uint32_t t0 = timer_hw->timerawl;

  // ---- 1. Latch what was computed last tick.  Fixed-length work => constant UD timing. ----
  if (g_ftwPending) { writeFTWWords(g_ftwHi, g_ftwLo); g_ftwPending = false; }
  writePhaseAmp(s_nextPh, s_nextAmp);
  pulseUpdate();

  // ---- 2. Next input sample (every INTERP ticks) + rate-lock bookkeeping ----
  if (s_sub == 0) {
    s_z0i = s_z1i; s_z0q = s_z1q;
    uint16_t head = s_head, tail = s_tail;
    int fill = (head - tail) & RING_MASK;
    s_fillAvg += ((float)fill - s_fillAvg) * (1.0f / 1024.0f);
    g_fillAvg = s_fillAvg;

    if (!s_running && fill >= RING_TARGET) s_running = true;
    if (s_running) {
      if (fill > 0) {
        uint32_t w = s_ring[tail];
        __dmb();
        s_tail = (uint16_t)((tail + 1) & RING_MASK);
        bool ptt; decodeWord(w, s_z1i, s_z1q, ptt);
        s_pttBit = ptt;
        g_lastI = (int16_t)s_z1i; g_lastQ = (int16_t)s_z1q; g_lastPtt = ptt;
      } else {
        g_underruns++;          // hold last sample; re-prefill before resuming
        s_running = false;
      }
    }

    // Steer the update period: FIFO fuller than target -> shorter period (run faster).
    int32_t nomQ16 = (int32_t)(s_periodNomCounts << 16);
    if (g_testActive) {
      s_periodQ16 = nomQ16;
    } else {
      float err = s_fillAvg - (float)RING_TARGET;
      int32_t adj = (int32_t)(err * PLL_GAIN_Q16);
      int32_t lim = nomQ16 / 100;                       // +-1 %
      if (adj > lim) adj = lim;
      if (adj < -lim) adj = -lim;
      s_periodQ16 = nomQ16 - adj;
    }
  }

  // ---- 3. Compute the registers for the NEXT tick ----
  bool localPtt = (gpio_get(PIN_LOCAL_PTT) == 0);
  float zi, zq;
  bool txOk;
  if (g_testActive) {
    s_testPh += 6.2831853f * g_testHz / (float)(SAMPLE_RATE_HZ * INTERP);
    if (s_testPh > 6.2831853f) s_testPh -= 6.2831853f;
    if (s_testPh < 0.0f) s_testPh += 6.2831853f;
    zi = (float)g_testAmp * cosf(s_testPh);
    zq = (float)g_testAmp * sinf(s_testPh);
    txOk = localPtt;
  } else {
    float f = (float)s_sub * (1.0f / (float)INTERP);
    zi = s_z0i + (s_z1i - s_z0i) * f;        // linear interpolation of I/Q (not of polar values!)
    zq = s_z0q + (s_z1q - s_z0q) * f;
    bool fresh = (time_us_32() - g_lastPacketUs) < STALE_TIMEOUT_US;
    txOk = localPtt && fresh && s_pttBit && s_running;
  }

  float A = sqrtf(zi * zi + zq * zq);
  if (!txOk) {
    s_nextAmp = 0;                              // carrier parked, phase unchanged
    s_nextPh = s_lastPh;
  } else {
    if (A > 4095.0f) A = 4095.0f;
    s_nextAmp = (uint16_t)(A + 0.5f);
    float ph = fastAtan2(zq, zi);               // -pi..pi
    float v = ph * (16384.0f / 6.2831853f);
    int32_t w = (int32_t)(v + (v >= 0.0f ? 0.5f : -0.5f));
    s_lastPh = (uint16_t)(w & 0x3FFF);
    s_nextPh = s_lastPh;
  }

  if (++s_sub >= INTERP) s_sub = 0;

  // ---- 4. Program the next period (takes effect at the coming wrap), 1-count dithering ----
  s_acc += (uint32_t)s_periodQ16;                   // Q16 clk_sys counts per DDS update
  uint32_t counts = s_acc >> 16;
  s_acc &= 0xFFFF;
  pwm_set_wrap(PWM_SLICE, counts - 1);
  g_periodCounts = counts;

  uint32_t dt = timer_hw->timerawl - t0;
  if (dt > g_isrMaxUs) g_isrMaxUs = dt;
}

// ---- core1: UART1 receiver, self-synchronising on the MSB-set byte ----
void core1_entry() {
  gpio_set_function(PIN_UART_TX, GPIO_FUNC_UART);
  gpio_set_function(PIN_UART_RX, GPIO_FUNC_UART);
  uart_init(uart1, LINK_BAUD);
  uart_set_format(uart1, 8, 1, UART_PARITY_NONE);
  uart_set_fifo_enabled(uart1, true);

  uint32_t w = 0;
  int idx = 0;
  while (true) {
    if (!uart_is_readable(uart1)) continue;
    uint32_t dr = uart_get_hw(uart1)->dr;
    uint8_t b = dr & 0xFF;
    if (dr & 0xF00) { idx = 0; g_badPackets++; continue; }      // framing/parity/overrun/break
    if (b & 0x80) {
      w = b & 0x7F; idx = 1;                                    // start of a packet
    } else if (idx > 0) {
      w = (w << 7) | b;
      if (++idx == 4) {
        idx = 0;
        if (__builtin_popcount(w) & 1) { g_badPackets++; continue; }
        uint16_t head = s_head;
        uint16_t next = (uint16_t)((head + 1) & RING_MASK);
        if (next == s_tail) { g_ringFull++; continue; }
        s_ring[head] = w;
        __dmb();
        s_head = next;
        g_lastPacketUs = time_us_32();
        g_packetCount++;
      }
    }
  }
}

void setup() {
  Serial.begin(115200);

  for (int p = 0; p <= 17; p++) { gpio_init(p); gpio_set_dir(p, GPIO_OUT); gpio_put(p, 0); }
  gpio_put(PIN_WR, 1); gpio_put(PIN_RD, 1); gpio_put(PIN_SP, 1);   // parallel mode
  gpio_init(PIN_UD); gpio_set_dir(PIN_UD, GPIO_OUT); gpio_put(PIN_UD, 0);
  gpio_init(PIN_LOCAL_PTT); gpio_set_dir(PIN_LOCAL_PTT, GPIO_IN); gpio_pull_up(PIN_LOCAL_PTT);

  gpio_put(PIN_RST, 1); delay(10); gpio_put(PIN_RST, 0); delay(10);

  writeReg(0x1D, 0x10);            // comparator off, both DACs on
  writeReg(0x1E, 0x40 | MULT);     // PLL on, x10 -> 300 MHz
  writeReg(0x20, 0x20);            // OSK enabled, manual multiplier
  writeReg(0x1F, 0x00);            // external update clock (UD pin), single-tone mode
  uint32_t hi, lo;
  computeFTW(g_carrierHz, hi, lo);
  writeFTWWords(hi, lo);
  writePhaseAmp(0, 0);
  pulseUpdate();

  s_periodNomCounts = (uint32_t)((double)clock_get_hz(clk_sys) / (SAMPLE_RATE_HZ * INTERP) + 0.5);
  s_periodNomUs = (uint32_t)(1e6 / (SAMPLE_RATE_HZ * INTERP)) + 1;
  s_periodQ16 = (int32_t)(s_periodNomCounts << 16);

  delay(1000);
  Serial.println("AD9854 SSB driver v3 up (phase-register polar modulation). Starts MUTED.");
  Serial.printf("DDS update rate %.1f kHz (INTERP=%d), clk_sys %.1f MHz, %lu counts/update\n",
                SAMPLE_RATE_HZ * INTERP / 1e3, INTERP, clock_get_hz(clk_sys) / 1e6, (unsigned long)s_periodNomCounts);
  Serial.println("Commands: <MHz> carrier, x<Hz>[,amp] single-tone test (x = off), i status");

  multicore_launch_core1(core1_entry);

  pwm_config c = pwm_get_default_config();
  pwm_config_set_clkdiv(&c, 1.0f);
  pwm_config_set_wrap(&c, s_periodNomCounts - 1);
  pwm_init(PWM_SLICE, &c, false);
  pwm_clear_irq(PWM_SLICE);
  pwm_set_irq0_enabled(PWM_SLICE, true);
  irq_set_exclusive_handler(PWM_DEFAULT_IRQ_NUM(), ddsTick);
  irq_set_priority(PWM_DEFAULT_IRQ_NUM(), 0);       // highest: beats USB for timing
  irq_set_enabled(PWM_DEFAULT_IRQ_NUM(), true);
  pwm_set_enabled(PWM_SLICE, true);
}

void loop() {
  if (!Serial.available()) return;
  char c = Serial.peek();
  if (c == 'i') {
    Serial.read();
    uint32_t age = time_us_32() - g_lastPacketUs;
    Serial.printf("carrier=%.6f MHz  packets=%lu  bad=%lu  ringFull=%lu  underruns=%lu\n",
                  g_carrierHz / 1e6, (unsigned long)g_packetCount, (unsigned long)g_badPackets,
                  (unsigned long)g_ringFull, (unsigned long)g_underruns);
    Serial.printf("fifo fill avg=%.1f (target %d)  update period=%lu counts (nominal %lu)\n",
                  g_fillAvg, RING_TARGET, (unsigned long)g_periodCounts, (unsigned long)s_periodNomCounts);
    Serial.printf("isr max=%lu us (budget %lu us)  last I=%d Q=%d  packetPTT=%s localPTT=%s  %s\n",
                  (unsigned long)g_isrMaxUs, (unsigned long)s_periodNomUs, g_lastI, g_lastQ,
                  g_lastPtt ? "TX" : "off", gpio_get(PIN_LOCAL_PTT) == 0 ? "TX" : "off",
                  age < STALE_TIMEOUT_US ? "fresh" : "STALE");
    Serial.printf("test tone=%s (%.1f Hz, amp %u)\n", g_testActive ? "ON" : "off", g_testHz, g_testAmp);
    g_isrMaxUs = 0;
  } else if (c == 'x') {
    Serial.read();
    int p = Serial.peek();
    if (!isdigit(p) && p != '-' && p != '.') {
      g_testActive = false;
      Serial.println("Test tone OFF -- back to Teensy data.");
    } else {
      g_testHz = Serial.parseFloat();
      if (Serial.peek() == ',') {
        Serial.read();
        int a = Serial.parseInt();
        if (a < 0) a = 0;
        if (a > 4095) a = 4095;
        g_testAmp = (uint16_t)a;
      }
      g_testActive = true;
      Serial.printf("Test tone ON: carrier%+.1f Hz, amp %u (local PTT required)\n", g_testHz, g_testAmp);
    }
  } else if (isdigit(c) || c == '.') {
    float mhz = Serial.parseFloat();
    if (mhz > 0.0f && mhz <= 120.0f) {
      g_carrierHz = mhz * 1e6;
      uint32_t hi, lo;
      computeFTW(g_carrierHz, hi, lo);
      g_ftwHi = hi; g_ftwLo = lo;
      __dmb();
      g_ftwPending = true;                          // the ISR does the bus writes, never loop()
      Serial.printf("Carrier = %.6f MHz\n", mhz);
    }
  } else {
    Serial.read();
  }
}
