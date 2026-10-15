// ============================================================
//  AudioEffectCompressor.h
//  A single-channel feed-forward RMS compressor for the
//  Teensy Audio Library (Teensyduino).
//
//  Drop this file in the same folder as your .ino sketch.
//  No other changes to the sketch are needed.
//
//  API:
//    compressor.threshold(dBFS);   // e.g. -18.0f
//    compressor.ratio(N);          // e.g. 4.0f  (N:1)
//    compressor.attack(ms);        // e.g. 5.0f
//    compressor.release(ms);       // e.g. 150.0f
// ============================================================

#pragma once
#include <Arduino.h>
#include <AudioStream.h>
#include <arm_math.h>   // CMSIS DSP — available on all Teensy 3.x / 4.x

class AudioEffectCompressor : public AudioStream {
public:
    AudioEffectCompressor()
        : AudioStream(1, inputQueueArray),
          thresholdLinear(0.1259f),   // -18 dBFS default
          slope(0.75f),               // ratio 4:1  → slope = 1 - 1/R
          attackCoef(0.0f),
          releaseCoef(0.0f),
          envelope(0.0f)
    {
        // Set sensible defaults so the object works even if the
        // user never calls the setters.
        setAttack(5.0f);
        setRelease(150.0f);
    }

    // ── Setter API (match the sketch's call-style) ────────────
    void threshold(float dBFS) {
        // Convert dBFS → linear amplitude  (|dBFS| / 20)
        thresholdLinear = powf(10.0f, dBFS / 20.0f);
    }

    void ratio(float r) {
        if (r < 1.0f) r = 1.0f;
        // slope used in gain-computer:  outputDB = threshold + (inputDB - threshold) / R
        // expressed in linear gain as   (1 - 1/R)
        slope = 1.0f - (1.0f / r);
    }

    void attack(float ms)  { setAttack(ms); }
    void release(float ms) { setRelease(ms); }

    // ── AudioStream interface ─────────────────────────────────
    void update() override {
        audio_block_t *block = receiveWritable(0);
        if (!block) return;

        //const float fs = AUDIO_SAMPLE_RATE_EXACT;
        int16_t *p = block->data;

        for (int i = 0; i < AUDIO_BLOCK_SAMPLES; i++) {
            // 1. Convert sample to float [-1, +1]
            float x = (float)p[i] / 32768.0f;

            // 2. Level detection: full-wave rectifier → 1-pole envelope follower
            float absX = fabsf(x);
            float coef = (absX > envelope) ? attackCoef : releaseCoef;
            envelope = absX + coef * (envelope - absX);

            // 3. Gain computer (in dB domain)
            float gainDB = 0.0f;
            if (envelope > thresholdLinear) {
                float overDB = 20.0f * log10f(envelope / thresholdLinear);
                gainDB = -slope * overDB;   // negative → attenuation
            }

            // 4. Apply gain
            float gainLinear = powf(10.0f, gainDB / 20.0f);
            float y = x * gainLinear;

            // 5. Soft clip & convert back to int16
            if      (y >  1.0f) y =  1.0f;
            else if (y < -1.0f) y = -1.0f;
            p[i] = (int16_t)(y * 32767.0f);
        }

        transmit(block);
        AudioStream::release(block);
    }

private:
    audio_block_t *inputQueueArray[1];

    float thresholdLinear;
    float slope;
    float attackCoef;
    float releaseCoef;
    float envelope;

    // τ = -1 / (fs * ln(1 - coef))  ↔  coef = 1 - e^(-1/(fs*τ_samples))
    void setAttack(float ms) {
        float tau = ms * 0.001f * AUDIO_SAMPLE_RATE_EXACT;
        attackCoef  = (tau > 0.0f) ? expf(-1.0f / tau) : 0.0f;
    }
    void setRelease(float ms) {
        float tau = ms * 0.001f * AUDIO_SAMPLE_RATE_EXACT;
        releaseCoef = (tau > 0.0f) ? expf(-1.0f / tau) : 0.0f;
    }
};
