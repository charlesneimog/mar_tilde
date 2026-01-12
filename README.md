# `mar~`

`mar~` is a Pure Data external for loading and playing audio files with automatic resampling.

## Features

- Pure Data external (`mar~`)
- Audio format support:
  - MP3 (via minimp3)
  - WAV and AIFF (via AudioFile)
  - FLAC (TODO)
- High-quality resampling (via r8brain)

## Building

```sh
cmake . -B build
cmake --build build
```

## Dependencies

* [minimp3](https://github.com/lieff/minimp3)
* [AudioFile](https://github.com/adamstark/AudioFile)
* [r8brain-free-src](https://github.com/avaneev/r8brain-free-src)

## Status

* FLAC support: TODO
