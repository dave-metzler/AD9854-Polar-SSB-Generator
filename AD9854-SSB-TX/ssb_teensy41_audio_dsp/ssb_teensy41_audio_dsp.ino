// SSB polar modulator -- AUDIO DSP STAGE v3, for Teensy 4.1 + Teensy Audio Shield (SGTL5000).
// UNTESTED on hardware (the back-end design was verified in a Python model).
//
// Part 1 of 2: capture mic/line audio, compress, band-limit, run the +-45 deg dual-Hilbert FIR pair,
// recover a standard (I,Q) pair and stream it, sample by sample, over a 4 Mbaud UART to the Pico 2
// (ssb_pico2_dds_driver.ino).  The Pico converts to polar and drives the AD9854 phase/amplitude
// registers.
//
// REVISION 9 (this rewrite) -- why the old polar-on-Teensy scheme was dropped
// ---------------------------------------------------------------------------
// Modelled end to end with speech-like input, the v2 chain gave ~25 dB opposite-sideband rejection and
// a carrier ~30 dB below the wanted sideband, while this file's Hilbert/bandpass front end alone is
// good for ~60 dB.  The loss was in the polar back end, which is a sampled system:
//   * envelope/frequency were decimated 5:1 to 8.82 kHz (sampling images ~15 dB down, 5.8 kHz away);
//   * the DDS was driven by FREQUENCY, i.e. it integrates -- so the envelope gate, hysteresis, one-pole
//     smoother, envelope-weighted decimation average, +-16 kHz clamp and the Pico's 100 Hz/cycle slew
//     limiter (all added to fight envelope-null glitches) each corrupted the integrated phase for good;
//   * packet pacing used "sincePacketSend = 0", which makes the average send rate slightly SLOWER than
//     8820 Hz, so the 256-entry ring filled up and dropped samples periodically.
// Now this stage does only linear, well-behaved math: I/Q at the audio rate, no gating, no decimation.
// The Pico takes I/Q to polar and sets the DDS phase register directly (no integration), updating at
// 88.2 kHz with the update rate locked to this board's audio clock via FIFO fill level.
//
// Also new: an explicit band limit BEFORE the Hilbert pair (8th-order Butterworth LP at 3 kHz + 8th-order
// HP at 250 Hz).  The Q=0.4 post-Hilbert bandpass alone is only ~20 dB down at 6 kHz (a real SSB signal
// needs far better) and does nothing about sub-300 Hz content, where the 200-tap Hilbert pair's image
// rejection falls to 22 dB at 300 Hz and ~7 dB at 100 Hz.  The Q=0.4 post-filter is kept as the voice
// tilt/shaping you tuned; it is no longer asked to be the band limiter.
//
// CREDIT: compressor and the 200-tap +-45 deg Hilbert tables are the user's own field-tested blocks.
//
// ---- Wiring ----
//   Teensy Audio Shield stacks on the Teensy 4.1.
//   PTT switch: one side GND, other to PTT_PIN (pullup, LOW = transmit); second pole to the Pico 2's
//     local PTT pin (GP22) as a redundant hardware mute.
//   UART to Pico 2: TX1 (pin 1) -> Pico GP21 (UART1 RX).  RX1 (pin 0) <- GP20 (unused).  Common GND.
//     Keep the link wire short (< 15 cm) -- it runs at 4 Mbaud.  If you see "bad" packets on the Pico
//     status, drop LINK_BAUD to 3000000 on BOTH boards.
//
// Serial commands (115200 baud, USB):
//   g<dB>   input gain trim, pre-compressor            t<dB> comp threshold   c<N> ratio   a<ms> attack
//   e<ms>   comp release                               k<x> comp makeup gain  b0/b1 comp active/bypass
//   n<dB>   SGTL5000 MIC preamp gain, 0-63 dB (default 20) -- the control to use with src1
//   r<0-15> SGTL5000 line-in level (default 5; LOWER = less sensitive; use if 'i' shows ADC peak ~1.0)
//   m<x>    post-filter gain into I/Q (default 8): raise until 'i' shows peak |z| ~0.6-0.9, never >1
//   w<Hz>   pre-Hilbert low-pass corner (default 3000) d<Hz> pre-Hilbert high-pass corner (default 250)
//   o<Hz>   post-Hilbert bandpass centre (1500)        q<x> its Q (0.4)
//   src0/src1 line / mic input      u / l  USB / LSB (default USB)
//   h0/h1/h2 headphone monitor off/RAW/PROCESSED       hv<x> headphone volume
//   i       status

#include <Audio.h>
#include <Wire.h>
#include <math.h>
#include "AudioEffectCompressor.h"

// ======================= Tunable parameters =======================
#define NO_HILBERT_COEFFS 200
const int PTT_PIN = 2;                // Teensy digital pin, internal pullup, active low
const uint32_t LINK_BAUD = 4000000;   // must match the Pico 2 side

// ======================= ±45 deg Hilbert FIR coefficients =======================
// Field-validated: N=200, Kaiser beta=8.6 (~80 dB stopband), fs=44100 Hz, wideband.
// Verified 90.000 deg differential phase from 500 Hz - 22 kHz. The upstream bandpass
// restricts the signal to 300-3000 Hz before it reaches these filters.
const short Hilbert_Plus_45_Coeffs[NO_HILBERT_COEFFS] = {
  (short)(32768 * -0.00000636189699459289),
  (short)(32768 * +0.00000000000000000018),
  (short)(32768 * -0.00001222392129728388),
  (short)(32768 * +0.00000000000000000014),
  (short)(32768 * -0.00002042495474875781),
  (short)(32768 * -0.00000000000000000006),
  (short)(32768 * -0.00003151175820731277),
  (short)(32768 * +0.00000000000000000060),
  (short)(32768 * -0.00004610689635467481),
  (short)(32768 * +0.00000000000000000025),
  (short)(32768 * -0.00006491337235708785),
  (short)(32768 * -0.00000000000000000048),
  (short)(32768 * -0.00008871905447125542),
  (short)(32768 * +0.00000000000000000119),
  (short)(32768 * -0.00011840088470155020),
  (short)(32768 * +0.00000000000000000009),
  (short)(32768 * -0.00015492887648053781),
  (short)(32768 * +0.00000000000000000323),
  (short)(32768 * -0.00019936992906629533),
  (short)(32768 * +0.00000000000000000170),
  (short)(32768 * -0.00025289151147945193),
  (short)(32768 * +0.00000000000000000314),
  (short)(32768 * -0.00031676529903824186),
  (short)(32768 * +0.00000000000000000008),
  (short)(32768 * -0.00039237088179112873),
  (short)(32768 * +0.00000000000000000166),
  (short)(32768 * -0.00048119970757343851),
  (short)(32768 * +0.00000000000000000380),
  (short)(32768 * -0.00058485947458829983),
  (short)(32768 * +0.00000000000000000688),
  (short)(32768 * -0.00070507925142179512),
  (short)(32768 * -0.00000000000000000032),
  (short)(32768 * -0.00084371567906035827),
  (short)(32768 * +0.00000000000000000291),
  (short)(32768 * -0.00100276070359690199),
  (short)(32768 * +0.00000000000000000720),
  (short)(32768 * -0.00118435140509196577),
  (short)(32768 * +0.00000000000000001294),
  (short)(32768 * -0.00139078263463435301),
  (short)(32768 * -0.00000000000000000097),
  (short)(32768 * -0.00162452335788858151),
  (short)(32768 * +0.00000000000000000485),
  (short)(32768 * -0.00188823784305856574),
  (short)(32768 * +0.00000000000000001230),
  (short)(32768 * -0.00218481314352544703),
  (short)(32768 * +0.00000000000000002233),
  (short)(32768 * -0.00251739473775600131),
  (short)(32768 * -0.00000000000000000324),
  (short)(32768 * -0.00288943274054190942),
  (short)(32768 * +0.00000000000000000647),
  (short)(32768 * -0.00330474184682656063),
  (short)(32768 * +0.00000000000000001974),
  (short)(32768 * -0.00376757919525306552),
  (short)(32768 * +0.00000000000000003624),
  (short)(32768 * -0.00428274576627975245),
  (short)(32768 * -0.00000000000000000777),
  (short)(32768 * -0.00485571894481148184),
  (short)(32768 * +0.00000000000000000841),
  (short)(32768 * -0.00549282676378390431),
  (short)(32768 * +0.00000000000000003042),
  (short)(32768 * -0.00620147854528216129),
  (short)(32768 * +0.00000000000000000971),
  (short)(32768 * -0.00699047287135014239),
  (short)(32768 * +0.00000000000000003689),
  (short)(32768 * -0.00787041318438902980),
  (short)(32768 * +0.00000000000000000971),
  (short)(32768 * -0.00885427572019564617),
  (short)(32768 * +0.00000000000000004401),
  (short)(32768 * -0.00995819711378333021),
  (short)(32768 * +0.00000000000000001165),
  (short)(32768 * -0.01120258547297215378),
  (short)(32768 * +0.00000000000000005437),
  (short)(32768 * -0.01261371903810677882),
  (short)(32768 * +0.00000000000000000906),
  (short)(32768 * -0.01422609946008015352),
  (short)(32768 * +0.00000000000000006343),
  (short)(32768 * -0.01608600848173297182),
  (short)(32768 * +0.00000000000000000777),
  (short)(32768 * -0.01825705073971429365),
  (short)(32768 * +0.00000000000000007637),
  (short)(32768 * -0.02082910759736848055),
  (short)(32768 * +0.00000000000000000518),
  (short)(32768 * -0.02393342967401571042),
  (short)(32768 * +0.00000000000000000259),
  (short)(32768 * -0.02776941138998956921),
  (short)(32768 * -0.00000000000000000259),
  (short)(32768 * -0.03265516073248370815),
  (short)(32768 * -0.00000000000000000518),
  (short)(32768 * -0.03913083066298363166),
  (short)(32768 * -0.00000000000000001294),
  (short)(32768 * -0.04819246424769733222),
  (short)(32768 * -0.00000000000000002071),
  (short)(32768 * -0.06189997688219894822),
  (short)(32768 * -0.00000000000000004142),
  (short)(32768 * -0.08531162182054254084),
  (short)(32768 * +0.00000000000000003107),
  (short)(32768 * -0.13505086845193309264),
  (short)(32768 * +0.00000000000000002071),
  (short)(32768 * -0.31640820495547311575),
  (short)(32768 * +0.00000000000000008285),
  (short)(32768 * +0.95000000000000006661),
  (short)(32768 * +0.00000000000000002071),
  (short)(32768 * +0.18953510473353563093),
  (short)(32768 * +0.00000000000000004142),
  (short)(32768 * +0.10469681078582140921),
  (short)(32768 * +0.00000000000000009838),
  (short)(32768 * +0.07183338432886621217),
  (short)(32768 * +0.00000000000000007767),
  (short)(32768 * +0.05426095682104956769),
  (short)(32768 * +0.00000000000000007249),
  (short)(32768 * +0.04324646636649577547),
  (short)(32768 * +0.00000000000000006213),
  (short)(32768 * +0.03564696277933902119),
  (short)(32768 * +0.00000000000000005437),
  (short)(32768 * +0.03005430517890990763),
  (short)(32768 * +0.00000000000000004660),
  (short)(32768 * +0.02574370767407021393),
  (short)(32768 * +0.00000000000000004660),
  (short)(32768 * +0.02230430291894012693),
  (short)(32768 * +0.00000000000000004401),
  (short)(32768 * +0.01948595601987595724),
  (short)(32768 * +0.00000000000000004013),
  (short)(32768 * +0.01712777380016347953),
  (short)(32768 * -0.00000000000000001812),
  (short)(32768 * +0.01512162712710481138),
  (short)(32768 * +0.00000000000000003366),
  (short)(32768 * +0.01339218957593901228),
  (short)(32768 * -0.00000000000000001294),
  (short)(32768 * +0.01188537985005565754),
  (short)(32768 * +0.00000000000000003042),
  (short)(32768 * +0.01056135116246077407),
  (short)(32768 * -0.00000000000000000712),
  (short)(32768 * +0.00939006950760831044),
  (short)(32768 * +0.00000000000000002589),
  (short)(32768 * +0.00834843082814038309),
  (short)(32768 * -0.00000000000000000453),
  (short)(32768 * +0.00741832732937679901),
  (short)(32768 * +0.00000000000000002265),
  (short)(32768 * +0.00658531825359155634),
  (short)(32768 * -0.00000000000000000194),
  (short)(32768 * +0.00583769656306775634),
  (short)(32768 * -0.00000000000000002039),
  (short)(32768 * +0.00516582145392418396),
  (short)(32768 * +0.00000000000000003430),
  (short)(32768 * +0.00456163334893700586),
  (short)(32768 * +0.00000000000000001586),
  (short)(32768 * +0.00401829665243560921),
  (short)(32768 * +0.00000000000000000097),
  (short)(32768 * +0.00352993355506319275),
  (short)(32768 * -0.00000000000000001068),
  (short)(32768 * +0.00309142376160698390),
  (short)(32768 * +0.00000000000000002201),
  (short)(32768 * +0.00269825262753427058),
  (short)(32768 * +0.00000000000000001019),
  (short)(32768 * +0.00234639528812497248),
  (short)(32768 * +0.00000000000000000146),
  (short)(32768 * +0.00203222783893897661),
  (short)(32768 * -0.00000000000000000518),
  (short)(32768 * +0.00175245903353380209),
  (short)(32768 * +0.00000000000000001327),
  (short)(32768 * +0.00150407765747069999),
  (short)(32768 * +0.00000000000000000639),
  (short)(32768 * +0.00128431194579424804),
  (short)(32768 * +0.00000000000000000154),
  (short)(32768 * +0.00109059828522910479),
  (short)(32768 * -0.00000000000000000202),
  (short)(32768 * +0.00092055708319910691),
  (short)(32768 * +0.00000000000000000728),
  (short)(32768 * +0.00077197416185649308),
  (short)(32768 * +0.00000000000000000368),
  (short)(32768 * +0.00064278639359589861),
  (short)(32768 * +0.00000000000000000097),
  (short)(32768 * +0.00053107056767095313),
  (short)(32768 * -0.00000000000000000077),
  (short)(32768 * +0.00043503468845665735),
  (short)(32768 * +0.00000000000000000352),
  (short)(32768 * +0.00035301107086852771),
  (short)(32768 * +0.00000000000000000176),
  (short)(32768 * +0.00028345072908744919),
  (short)(32768 * +0.00000000000000000059),
  (short)(32768 * +0.00022491865943263372),
  (short)(32768 * +0.00000000000000000259),
  (short)(32768 * +0.00017608970306065260),
  (short)(32768 * -0.00000000000000000071),
  (short)(32768 * +0.00013574474359599197),
  (short)(32768 * +0.00000000000000000074),
  (short)(32768 * +0.00010276705208927049),
  (short)(32768 * +0.00000000000000000150),
  (short)(32768 * +0.00007613863927559657),
  (short)(32768 * -0.00000000000000000005),
  (short)(32768 * +0.00005493651477085162),
  (short)(32768 * +0.00000000000000000047),
  (short)(32768 * +0.00003832878595441639),
  (short)(32768 * -0.00000000000000000024),
  (short)(32768 * +0.00002557055687064513),
  (short)(32768 * +0.00000000000000000007),
  (short)(32768 * +0.00001599961032998473),
  (short)(32768 * +0.00000000000000000017),
  (short)(32768 * +0.00000903187512250060),
  (short)(32768 * -0.00000000000000000002)
};

const short Hilbert_Minus_45_Coeffs[NO_HILBERT_COEFFS] = {
  (short)(32768 * -0.00000000000000000002),
  (short)(32768 * +0.00000903187512250060),
  (short)(32768 * +0.00000000000000000017),
  (short)(32768 * +0.00001599961032998473),
  (short)(32768 * +0.00000000000000000007),
  (short)(32768 * +0.00002557055687064513),
  (short)(32768 * -0.00000000000000000024),
  (short)(32768 * +0.00003832878595441639),
  (short)(32768 * +0.00000000000000000047),
  (short)(32768 * +0.00005493651477085162),
  (short)(32768 * -0.00000000000000000005),
  (short)(32768 * +0.00007613863927559657),
  (short)(32768 * +0.00000000000000000150),
  (short)(32768 * +0.00010276705208927049),
  (short)(32768 * +0.00000000000000000074),
  (short)(32768 * +0.00013574474359599197),
  (short)(32768 * -0.00000000000000000071),
  (short)(32768 * +0.00017608970306065260),
  (short)(32768 * +0.00000000000000000259),
  (short)(32768 * +0.00022491865943263372),
  (short)(32768 * +0.00000000000000000059),
  (short)(32768 * +0.00028345072908744919),
  (short)(32768 * +0.00000000000000000176),
  (short)(32768 * +0.00035301107086852771),
  (short)(32768 * +0.00000000000000000352),
  (short)(32768 * +0.00043503468845665735),
  (short)(32768 * -0.00000000000000000077),
  (short)(32768 * +0.00053107056767095313),
  (short)(32768 * +0.00000000000000000097),
  (short)(32768 * +0.00064278639359589861),
  (short)(32768 * +0.00000000000000000368),
  (short)(32768 * +0.00077197416185649308),
  (short)(32768 * +0.00000000000000000728),
  (short)(32768 * +0.00092055708319910691),
  (short)(32768 * -0.00000000000000000202),
  (short)(32768 * +0.00109059828522910479),
  (short)(32768 * +0.00000000000000000154),
  (short)(32768 * +0.00128431194579424804),
  (short)(32768 * +0.00000000000000000639),
  (short)(32768 * +0.00150407765747069999),
  (short)(32768 * +0.00000000000000001327),
  (short)(32768 * +0.00175245903353380209),
  (short)(32768 * -0.00000000000000000518),
  (short)(32768 * +0.00203222783893897661),
  (short)(32768 * +0.00000000000000000146),
  (short)(32768 * +0.00234639528812497248),
  (short)(32768 * +0.00000000000000001019),
  (short)(32768 * +0.00269825262753427058),
  (short)(32768 * +0.00000000000000002201),
  (short)(32768 * +0.00309142376160698390),
  (short)(32768 * -0.00000000000000001068),
  (short)(32768 * +0.00352993355506319275),
  (short)(32768 * +0.00000000000000000097),
  (short)(32768 * +0.00401829665243560921),
  (short)(32768 * +0.00000000000000001586),
  (short)(32768 * +0.00456163334893700586),
  (short)(32768 * +0.00000000000000003430),
  (short)(32768 * +0.00516582145392418396),
  (short)(32768 * -0.00000000000000002039),
  (short)(32768 * +0.00583769656306775634),
  (short)(32768 * -0.00000000000000000194),
  (short)(32768 * +0.00658531825359155634),
  (short)(32768 * +0.00000000000000002265),
  (short)(32768 * +0.00741832732937679901),
  (short)(32768 * -0.00000000000000000453),
  (short)(32768 * +0.00834843082814038309),
  (short)(32768 * +0.00000000000000002589),
  (short)(32768 * +0.00939006950760831044),
  (short)(32768 * -0.00000000000000000712),
  (short)(32768 * +0.01056135116246077407),
  (short)(32768 * +0.00000000000000003042),
  (short)(32768 * +0.01188537985005565754),
  (short)(32768 * -0.00000000000000001294),
  (short)(32768 * +0.01339218957593901228),
  (short)(32768 * +0.00000000000000003366),
  (short)(32768 * +0.01512162712710481138),
  (short)(32768 * -0.00000000000000001812),
  (short)(32768 * +0.01712777380016347953),
  (short)(32768 * +0.00000000000000004013),
  (short)(32768 * +0.01948595601987595724),
  (short)(32768 * +0.00000000000000004401),
  (short)(32768 * +0.02230430291894012693),
  (short)(32768 * +0.00000000000000004660),
  (short)(32768 * +0.02574370767407021393),
  (short)(32768 * +0.00000000000000004660),
  (short)(32768 * +0.03005430517890990763),
  (short)(32768 * +0.00000000000000005437),
  (short)(32768 * +0.03564696277933902119),
  (short)(32768 * +0.00000000000000006213),
  (short)(32768 * +0.04324646636649577547),
  (short)(32768 * +0.00000000000000007249),
  (short)(32768 * +0.05426095682104956769),
  (short)(32768 * +0.00000000000000007767),
  (short)(32768 * +0.07183338432886621217),
  (short)(32768 * +0.00000000000000009838),
  (short)(32768 * +0.10469681078582140921),
  (short)(32768 * +0.00000000000000004142),
  (short)(32768 * +0.18953510473353563093),
  (short)(32768 * +0.00000000000000002071),
  (short)(32768 * +0.95000000000000006661),
  (short)(32768 * +0.00000000000000008285),
  (short)(32768 * -0.31640820495547311575),
  (short)(32768 * +0.00000000000000002071),
  (short)(32768 * -0.13505086845193309264),
  (short)(32768 * +0.00000000000000003107),
  (short)(32768 * -0.08531162182054254084),
  (short)(32768 * -0.00000000000000004142),
  (short)(32768 * -0.06189997688219894822),
  (short)(32768 * -0.00000000000000002071),
  (short)(32768 * -0.04819246424769733222),
  (short)(32768 * -0.00000000000000001294),
  (short)(32768 * -0.03913083066298363166),
  (short)(32768 * -0.00000000000000000518),
  (short)(32768 * -0.03265516073248370815),
  (short)(32768 * -0.00000000000000000259),
  (short)(32768 * -0.02776941138998956921),
  (short)(32768 * +0.00000000000000000259),
  (short)(32768 * -0.02393342967401571042),
  (short)(32768 * +0.00000000000000000518),
  (short)(32768 * -0.02082910759736848055),
  (short)(32768 * +0.00000000000000007637),
  (short)(32768 * -0.01825705073971429365),
  (short)(32768 * +0.00000000000000000777),
  (short)(32768 * -0.01608600848173297182),
  (short)(32768 * +0.00000000000000006343),
  (short)(32768 * -0.01422609946008015352),
  (short)(32768 * +0.00000000000000000906),
  (short)(32768 * -0.01261371903810677882),
  (short)(32768 * +0.00000000000000005437),
  (short)(32768 * -0.01120258547297215378),
  (short)(32768 * +0.00000000000000001165),
  (short)(32768 * -0.00995819711378333021),
  (short)(32768 * +0.00000000000000004401),
  (short)(32768 * -0.00885427572019564617),
  (short)(32768 * +0.00000000000000000971),
  (short)(32768 * -0.00787041318438902980),
  (short)(32768 * +0.00000000000000003689),
  (short)(32768 * -0.00699047287135014239),
  (short)(32768 * +0.00000000000000000971),
  (short)(32768 * -0.00620147854528216129),
  (short)(32768 * +0.00000000000000003042),
  (short)(32768 * -0.00549282676378390431),
  (short)(32768 * +0.00000000000000000841),
  (short)(32768 * -0.00485571894481148184),
  (short)(32768 * -0.00000000000000000777),
  (short)(32768 * -0.00428274576627975245),
  (short)(32768 * +0.00000000000000003624),
  (short)(32768 * -0.00376757919525306552),
  (short)(32768 * +0.00000000000000001974),
  (short)(32768 * -0.00330474184682656063),
  (short)(32768 * +0.00000000000000000647),
  (short)(32768 * -0.00288943274054190942),
  (short)(32768 * -0.00000000000000000324),
  (short)(32768 * -0.00251739473775600131),
  (short)(32768 * +0.00000000000000002233),
  (short)(32768 * -0.00218481314352544703),
  (short)(32768 * +0.00000000000000001230),
  (short)(32768 * -0.00188823784305856574),
  (short)(32768 * +0.00000000000000000485),
  (short)(32768 * -0.00162452335788858151),
  (short)(32768 * -0.00000000000000000097),
  (short)(32768 * -0.00139078263463435301),
  (short)(32768 * +0.00000000000000001294),
  (short)(32768 * -0.00118435140509196577),
  (short)(32768 * +0.00000000000000000720),
  (short)(32768 * -0.00100276070359690199),
  (short)(32768 * +0.00000000000000000291),
  (short)(32768 * -0.00084371567906035827),
  (short)(32768 * -0.00000000000000000032),
  (short)(32768 * -0.00070507925142179512),
  (short)(32768 * +0.00000000000000000688),
  (short)(32768 * -0.00058485947458829983),
  (short)(32768 * +0.00000000000000000380),
  (short)(32768 * -0.00048119970757343851),
  (short)(32768 * +0.00000000000000000166),
  (short)(32768 * -0.00039237088179112873),
  (short)(32768 * +0.00000000000000000008),
  (short)(32768 * -0.00031676529903824186),
  (short)(32768 * +0.00000000000000000314),
  (short)(32768 * -0.00025289151147945193),
  (short)(32768 * +0.00000000000000000170),
  (short)(32768 * -0.00019936992906629533),
  (short)(32768 * +0.00000000000000000323),
  (short)(32768 * -0.00015492887648053781),
  (short)(32768 * +0.00000000000000000009),
  (short)(32768 * -0.00011840088470155020),
  (short)(32768 * +0.00000000000000000119),
  (short)(32768 * -0.00008871905447125542),
  (short)(32768 * -0.00000000000000000048),
  (short)(32768 * -0.00006491337235708785),
  (short)(32768 * +0.00000000000000000025),
  (short)(32768 * -0.00004610689635467481),
  (short)(32768 * +0.00000000000000000060),
  (short)(32768 * -0.00003151175820731277),
  (short)(32768 * -0.00000000000000000006),
  (short)(32768 * -0.00002042495474875781),
  (short)(32768 * +0.00000000000000000014),
  (short)(32768 * -0.00001222392129728388),
  (short)(32768 * +0.00000000000000000018),
  (short)(32768 * -0.00000636189699459289)
};

// ======================= I/Q extraction stage =======================
// p = plus45 branch, m = minus45 branch.  With p = A cos(PHI+45deg), m = A cos(PHI-45deg):
//   I = (p+m)/sqrt2 =  A cos(PHI),   Q = (p-m)/sqrt2 = -A sin(PHI)
// SIDEBAND SIGN (v3.1): bench testing on the air showed 'u' came out as LSB and 'l' as USB, i.e. the
// real p/m outputs are swapped relative to the tables' nominal labelling (the likely reason: CMSIS
// arm_fir expects its coefficient array time-reversed, and time-reversing the plus-45 table gives
// exactly the minus-45 table).  The fix is simply the sign of Q below.  'u' now means USB on the air.
// Linear math only: no atan2, no gating, no decimation.  Output is clamped by SCALING the vector
// (phase preserved) if |z| would exceed 1.0.
struct IQSample { int16_t i; int16_t q; };      // +-4095 == +-1.0 full scale

class SSBIQStage : public AudioStream {
 public:
  SSBIQStage() : AudioStream(2, inputQueueArray) {}
  virtual void update(void);

  void setSideband(int8_t sb) { sideband = sb; }       // +1 USB, -1 LSB

  bool pop(IQSample &out) {
    if (rbHead == rbTail) return false;
    out = ring[rbTail];
    rbTail = (rbTail + 1) & RB_MASK;
    return true;
  }
  uint16_t depth() const { return (rbHead - rbTail) & RB_MASK; }

  volatile uint32_t blockCount = 0;
  volatile uint32_t clipCount  = 0;     // samples where |z| > 1.0 had to be scaled down
  volatile uint32_t dropCount  = 0;     // ring full (Pico not draining) -- should stay 0
  volatile float peakMag = 0.0f;        // peak |z| since last resetPeak()
  volatile float lastMag = 0.0f;
  void resetPeak() { peakMag = 0.0f; }

 private:
  audio_block_t *inputQueueArray[2];   // [0]=plus45 (p), [1]=minus45 (m)
  int8_t sideband = +1;

  static const int RB_SIZE = 1024;
  static const int RB_MASK = RB_SIZE - 1;
  IQSample ring[RB_SIZE];
  volatile uint16_t rbHead = 0, rbTail = 0;
};

void SSBIQStage::update(void) {
  audio_block_t *pBlock = receiveReadOnly(0);
  audio_block_t *mBlock = receiveReadOnly(1);
  if (!pBlock || !mBlock) {
    if (pBlock) release(pBlock);
    if (mBlock) release(mBlock);
    blockCount++;
    return;
  }

  const float INV_SQRT2 = 0.70710678f;
  for (int n = 0; n < AUDIO_BLOCK_SAMPLES; n++) {
    float p = pBlock->data[n] / 32768.0f;
    float m = mBlock->data[n] / 32768.0f;
    float I = (p + m) * INV_SQRT2;
    float Q = -(float)sideband * (p - m) * INV_SQRT2;   // v3.1: sign flipped, see note above

    float mag = sqrtf(I * I + Q * Q);
    lastMag = mag;
    if (mag > peakMag) peakMag = mag;
    if (mag > 1.0f) { I /= mag; Q /= mag; clipCount++; }

    uint16_t next = (rbHead + 1) & RB_MASK;
    if (next == rbTail) { dropCount++; continue; }
    ring[rbHead].i = (int16_t)lroundf(I * 4095.0f);
    ring[rbHead].q = (int16_t)lroundf(Q * 4095.0f);
    rbHead = next;
  }

  blockCount++;
  release(pBlock);
  release(mBlock);
}

// ======================= Audio Library graph =======================
AudioInputI2S           micIn;
AudioAmplifier           inGainAmp;     // pre-compressor input trim ("g<dB>")
AudioEffectCompressor    comp;          // user-supplied, used as-is
AudioAmplifier           outGain;       // post-compressor makeup gain ("k<x>")
AudioFilterBiquad        preHP;         // NEW: 8th-order Butterworth high-pass (4 stages)
AudioFilterBiquad        preLP;         // NEW: 8th-order Butterworth low-pass  (4 stages)
AudioFilterFIR           hilbertPlus45;
AudioFilterFIR           hilbertMinus45;
AudioFilterBiquad        bandpassP;     // post-Hilbert matched bandpass, plus-45 branch
AudioFilterBiquad        bandpassM;     // post-Hilbert matched bandpass, minus-45 branch
AudioAmplifier           postGainP;
AudioAmplifier           postGainM;
SSBIQStage               ssbIQ;
AudioControlSGTL5000     sgtl5000;
AudioAnalyzePeak         pkAdc;         // level meters: any reading >= 0.98 means that stage is clipping
AudioAnalyzePeak         pkOut;
AudioAnalyzePeak         pkPre;
AudioMixer4              hpMix;
AudioOutputI2S           hpOut;

AudioConnection pc01(micIn,          0, inGainAmp,      0);
AudioConnection pc02(inGainAmp,      0, comp,           0);
AudioConnection pc03(comp,           0, outGain,        0);
AudioConnection pc03a(outGain,       0, preHP,          0);
AudioConnection pc03b(preHP,         0, preLP,          0);
AudioConnection pc04(preLP,          0, hilbertPlus45,  0);
AudioConnection pc05(preLP,          0, hilbertMinus45, 0);
AudioConnection pc06(hilbertPlus45,  0, bandpassP,      0);
AudioConnection pc07(hilbertMinus45, 0, bandpassM,      0);
AudioConnection pc08(bandpassP,      0, postGainP,      0);
AudioConnection pc09(bandpassM,      0, postGainM,      0);
AudioConnection pc10(postGainP,      0, ssbIQ,          0);
AudioConnection pc11(postGainM,      0, ssbIQ,          1);

AudioConnection pcHp1(inGainAmp,     0, hpMix,          0);   // RAW tap
AudioConnection pcHp2(postGainP,     0, hpMix,          1);   // PROCESSED tap
AudioConnection pcHp3(hpMix,         0, hpOut,          0);
AudioConnection pcHp4(hpMix,         0, hpOut,          1);

AudioConnection pcPk1(micIn,         0, pkAdc,          0);   // SGTL5000 ADC output
AudioConnection pcPk2(outGain,       0, pkOut,          0);   // after compressor + makeup gain (k)
AudioConnection pcPk3(preLP,         0, pkPre,          0);   // what actually enters the Hilbert FIRs

const int AUDIO_MEMORY_BLOCKS = 100;

bool pttActive = false;
bool compressorBypassed = false;
float savedThresholdDB = -28.0f;
float savedRatio = 4.0f;

// Post-Hilbert bandpass (user's proven shape).  Both branches ALWAYS get identical settings.
float bpCenterHz = 1500.0f;
float bpQ = 0.4f;

static void applyBandpass() {
  for (int s = 0; s < 4; s++) {
    bandpassP.setBandpass(s, bpCenterHz, bpQ);
    bandpassM.setBandpass(s, bpCenterHz, bpQ);
  }
}

// Pre-Hilbert band limit.  Butterworth 8th-order section Q values.
float preLpHz = 3000.0f;
float preHpHz = 250.0f;
static const float BW8_Q[4] = { 0.5098f, 0.6013f, 0.8999f, 2.5629f };

static void applyPreFilters() {
  for (int s = 0; s < 4; s++) {
    preLP.setLowpass(s, preLpHz, BW8_Q[s]);
    preHP.setHighpass(s, preHpHz, BW8_Q[s]);
  }
}

float postGainVal = 8.0f;
static void applyPostGain() {
  postGainP.gain(postGainVal);
  postGainM.gain(postGainVal);
}

// 4-byte, self-synchronising packet (see the Pico sketch header for the layout).
static inline void sendPacket(const IQSample &s, bool ptt) {
  uint32_t w = ((uint32_t)(ptt ? 1 : 0) << 27) |
               ((uint32_t)(s.i & 0x1FFF) << 13) |
               (uint32_t)(s.q & 0x1FFF);
  if (__builtin_popcount(w) & 1) w |= (1u << 26);            // even parity over the payload
  uint8_t b[4] = {
    (uint8_t)(0x80 | ((w >> 21) & 0x7F)),
    (uint8_t)((w >> 14) & 0x7F),
    (uint8_t)((w >> 7) & 0x7F),
    (uint8_t)(w & 0x7F)
  };
  Serial1.write(b, 4);
}

static uint8_t serial1TxBuf[1024];

void setup() {
  Serial.begin(115200);
  AudioMemory(AUDIO_MEMORY_BLOCKS);

  pinMode(PTT_PIN, INPUT_PULLUP);
  Serial1.addMemoryForWrite(serial1TxBuf, sizeof(serial1TxBuf));
  Serial1.begin(LINK_BAUD);

  sgtl5000.enable();
  sgtl5000.inputSelect(AUDIO_INPUT_LINEIN);   // default to line; "src1" switches to mic
  sgtl5000.lineInLevel(5);
  sgtl5000.micGain(20);
  sgtl5000.volume(0.5);

  hpMix.gain(0, 0.0f);
  hpMix.gain(1, 0.0f);

  inGainAmp.gain(1.0f);
  outGain.gain(4.0f);

  comp.threshold(savedThresholdDB);
  comp.ratio(savedRatio);
  comp.attack(5.0f);
  comp.release(150.0f);

  hilbertPlus45.begin (Hilbert_Plus_45_Coeffs,  NO_HILBERT_COEFFS);
  hilbertMinus45.begin(Hilbert_Minus_45_Coeffs, NO_HILBERT_COEFFS);

  applyPreFilters();
  applyBandpass();
  applyPostGain();

  delay(1500);
  Serial.println("SSB audio DSP v3 up: LINE IN, compressor active, USB, I/Q link to Pico at 4 Mbaud.");
  Serial.println("Commands: g t c a e k m b n r src0/src1 u/l w d o q h0/h1/h2 hv i");
  Serial.printf("Pre-Hilbert band limit: HP %.0f Hz, LP %.0f Hz (8th-order each)\n", preHpHz, preLpHz);
  Serial.printf("Post-filter gain m=%.2f -- raise until 'i' shows peak |z| ~0.6-0.9; clip must stay 0\n", postGainVal);
}

void loop() {
  pttActive = (digitalRead(PTT_PIN) == LOW);

  // Forward samples as fast as the UART takes them.  The Pico's FIFO + rate lock absorbs the 128-sample
  // bursts the audio library delivers; nothing here needs to be timed.
  IQSample s;
  while (Serial1.availableForWrite() >= 4 && ssbIQ.pop(s)) sendPacket(s, pttActive);

  if (Serial.available()) {
    char c = Serial.peek();
    if (c == 'i') {
      Serial.read();
      Serial.printf("blocks=%lu  |z| clipped=%lu  ringDrops=%lu  ring depth=%u/1024  PTT=%s  comp=%s\n",
                    ssbIQ.blockCount, ssbIQ.clipCount, ssbIQ.dropCount, ssbIQ.depth(),
                    pttActive ? "TX" : "off", compressorBypassed ? "bypassed" : "active");
      Serial.printf("peak |z| since last 'i' = %.3f  (target 0.6-0.9, must not exceed 1.0)  last=%.3f\n",
                    ssbIQ.peakMag, ssbIQ.lastMag);
      ssbIQ.resetPeak();
      Serial.printf("audio peaks since last 'i' (0-1, >=0.98 = CLIPPING): ADC=%.3f  post-comp/makeup=%.3f  into-Hilbert=%.3f\n",
                    pkAdc.available() ? pkAdc.read() : 0.0f,
                    pkOut.available() ? pkOut.read() : 0.0f,
                    pkPre.available() ? pkPre.read() : 0.0f);
      Serial.printf("CPU %.1f%% (peak %.1f%%)  mem %d/%d\n",
                    AudioProcessorUsage(), AudioProcessorUsageMax(),
                    AudioMemoryUsage(), AUDIO_MEMORY_BLOCKS);
      AudioProcessorUsageMaxReset();
    } else if (c == 'g') {
      Serial.read();
      float db = Serial.parseFloat();
      inGainAmp.gain(powf(10.0f, db / 20.0f));
      Serial.printf("Input gain = %.1f dB\n", db);
    } else if (c == 't') {
      Serial.read();
      savedThresholdDB = Serial.parseFloat();
      if (!compressorBypassed) comp.threshold(savedThresholdDB);
      Serial.printf("Threshold = %.1f dBFS\n", savedThresholdDB);
    } else if (c == 'c') {
      Serial.read();
      savedRatio = Serial.parseFloat();
      if (!compressorBypassed) comp.ratio(savedRatio);
      Serial.printf("Ratio = %.1f:1\n", savedRatio);
    } else if (c == 'a') {
      Serial.read();
      float ms = Serial.parseFloat();
      comp.attack(ms);
      Serial.printf("Attack = %.1f ms\n", ms);
    } else if (c == 'e') {
      Serial.read();
      float ms = Serial.parseFloat();
      comp.release(ms);
      Serial.printf("Release = %.1f ms\n", ms);
    } else if (c == 'k') {
      Serial.read();
      float g = Serial.parseFloat();
      outGain.gain(g);
      Serial.printf("Makeup gain = %.2fx\n", g);
    } else if (c == 'm') {
      Serial.read();
      postGainVal = Serial.parseFloat();
      applyPostGain();
      Serial.printf("Post-filter gain = %.2fx (watch 'i' peak |z|)\n", postGainVal);
    } else if (c == 'b') {
      Serial.read();
      int v = Serial.parseInt();
      compressorBypassed = (v != 0);
      if (compressorBypassed) { comp.threshold(0.0f); comp.ratio(1.0f); }
      else { comp.threshold(savedThresholdDB); comp.ratio(savedRatio); }
      Serial.printf("Compressor %s\n", compressorBypassed ? "BYPASSED" : "active");
    } else if (c == 'w') {
      Serial.read();
      preLpHz = Serial.parseFloat();
      applyPreFilters();
      Serial.printf("Pre-Hilbert LP = %.0f Hz\n", preLpHz);
    } else if (c == 'd') {
      Serial.read();
      preHpHz = Serial.parseFloat();
      applyPreFilters();
      Serial.printf("Pre-Hilbert HP = %.0f Hz\n", preHpHz);
    } else if (c == 'n') {
      Serial.read();
      int mg = Serial.parseInt();
      if (mg < 0) mg = 0;
      if (mg > 63) mg = 63;
      sgtl5000.micGain(mg);
      Serial.printf("Mic preamp gain = %d dB (applies to the MIC input; default 20). Lower it if 'i' shows ADC peak ~1.0\n", mg);
    } else if (c == 'r') {
      Serial.read();
      int lv = Serial.parseInt();
      if (lv < 0) lv = 0;
      if (lv > 15) lv = 15;
      sgtl5000.lineInLevel(lv);
      Serial.printf("Line-in level = %d  (0 = 3.12 Vpp full scale ... 5 = 1.33 Vpp ... 15 = 0.24 Vpp; LOWER number = less sensitive)\n", lv);
    } else if (c == 'o') {
      Serial.read();
      bpCenterHz = Serial.parseFloat();
      applyBandpass();
      Serial.printf("Bandpass center = %.0f Hz (Q=%.2f)\n", bpCenterHz, bpQ);
    } else if (c == 'q') {
      Serial.read();
      bpQ = Serial.parseFloat();
      applyBandpass();
      Serial.printf("Bandpass Q = %.2f (center=%.0f Hz)\n", bpQ, bpCenterHz);
    } else if (c == 'h') {
      Serial.read();
      char next = Serial.peek();
      if (next == 'v') {
        Serial.read();
        float v = Serial.parseFloat();
        if (v < 0.0f) v = 0.0f;
        if (v > 1.0f) v = 1.0f;
        sgtl5000.volume(v);
        Serial.printf("Headphone volume = %.2f\n", v);
      } else {
        int v = Serial.parseInt();
        if (v == 0) {
          hpMix.gain(0, 0.0f); hpMix.gain(1, 0.0f);
          Serial.println("Headphone monitor OFF");
        } else if (v == 1) {
          hpMix.gain(0, 1.0f); hpMix.gain(1, 0.0f);
          Serial.println("Headphone monitor = RAW (post-gain, pre-compressor)");
        } else if (v == 2) {
          hpMix.gain(0, 0.0f); hpMix.gain(1, 1.0f);
          Serial.println("Headphone monitor = PROCESSED (feeds the I/Q math)");
        }
      }
    } else if (c == 'u') {
      Serial.read();
      ssbIQ.setSideband(+1);
      Serial.println("Sideband = USB");
    } else if (c == 'l') {
      Serial.read();
      ssbIQ.setSideband(-1);
      Serial.println("Sideband = LSB");
    } else if (c == 's') {
      Serial.read();
      String tok = Serial.readStringUntil('\n');
      if (tok.startsWith("rc0")) { sgtl5000.inputSelect(AUDIO_INPUT_LINEIN); Serial.println("Input = LINE"); }
      else if (tok.startsWith("rc1")) { sgtl5000.inputSelect(AUDIO_INPUT_MIC); Serial.println("Input = MIC"); }
    } else {
      Serial.read();
    }
  }
}
