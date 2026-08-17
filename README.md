# `mar~`

`mar~` is a Pure Data external for loading and playing audio files with automatic resampling.

## Features

- Pure Data external (`mar~`)
- Audio format support:
  - MP3 (via minimp3)
  - WAV and AIFF (via AudioFile)
  - FLAC (TODO)
- High-quality resampling (via statically linked libsamplerate)

## Building

```sh
cmake . -B build
cmake --build build
```

## Dependencies

* [minimp3](https://github.com/lieff/minimp3)
* [AudioFile](https://github.com/adamstark/AudioFile)
* [libsamplerate](https://github.com/libsndfile/libsamplerate)

## Audio Example

* Example sound: [Piano loop: Sloka 3 - 90bpm 3/4](https://freesound.org/people/trader_one/sounds/683575/) by _trader_one_.
