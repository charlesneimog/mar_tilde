cmake_minimum_required(VERSION 3.25)
project(mar~)

# ╭──────────────────────────────────────╮
# │             CPM Package              │
# ╰──────────────────────────────────────╯
set(CPM_FILE ${CMAKE_BINARY_DIR}/CPM.cmake)
set(CPM_VERSION "0.42.0")
if(NOT EXISTS "${CPM_FILE}")
    file(DOWNLOAD "https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_VERSION}/CPM.cmake" ${CPM_FILE})
endif()
include(${CPM_FILE})

# ╭──────────────────────────────────────╮
# │               PD.CMAKE               │
# ╰──────────────────────────────────────╯
set(PDCMAKE_FILE ${CMAKE_BINARY_DIR}/pd.cmake)
set(PDCMAKE_VERSION "v0.2.11")
if(NOT EXISTS "${PDCMAKE_FILE}")
    file(DOWNLOAD
         # "https://github.com/pure-data/pd.cmake/releases/download/v2.0.11/pd.cmake"
         "https://raw.githubusercontent.com/charlesneimog/pd.cmake/refs/heads/flags-test/pd.cmake" ${PDCMAKE_FILE})
endif()
include(${PDCMAKE_FILE})

# ╭──────────────────────────────────────╮
# │                 MP3                  │
# ╰──────────────────────────────────────╯
cpmaddpackage("gh:lieff/minimp3#master")
if(minimp3_ADDED)
    include_directories("${minimp3_SOURCE_DIR}")
endif()

# ╭──────────────────────────────────────╮
# │              WAV & AIFF              │
# ╰──────────────────────────────────────╯
cpmaddpackage("gh:adamstark/AudioFile#master")

# ╭──────────────────────────────────────╮
# │                 FLAC                 │
# ╰──────────────────────────────────────╯
cpmaddpackage("gh:charlesneimog/miniflac#main")
if(miniflac_ADDED)
    include_directories("${miniflac_SOURCE_DIR}")
endif()

# ╭──────────────────────────────────────╮
# │              OGG Vorbis              │
# ╰──────────────────────────────────────╯
# stb_vorbis is a small, portable Ogg Vorbis decoder (not Opus).
cpmaddpackage("gh:nothings/stb#master")
if(stb_ADDED)
    include_directories("${stb_SOURCE_DIR}")
endif()

# ╭──────────────────────────────────────╮
# │              Resampler               │
# ╰──────────────────────────────────────╯
CPMAddPackage(
    NAME libsamplerate
    GITHUB_REPOSITORY libsndfile/libsamplerate
    GIT_TAG 0.2.2
    OPTIONS
        "BUILD_SHARED_LIBS OFF"
        "BUILD_TESTING OFF"
        "LIBSAMPLERATE_EXAMPLES OFF"
        "LIBSAMPLERATE_INSTALL OFF")

# The static library is linked into the Pd module, which requires PIC on ELF platforms.
set_target_properties(samplerate PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ╭──────────────────────────────────────╮
# │              Pd Object               │
# ╰──────────────────────────────────────╯
pd_add_external(
    mar~
    "${CMAKE_CURRENT_SOURCE_DIR}/mar~.cpp"
    LINK_LIBRARIES
    AudioFile
    SampleRate::samplerate)
pd_add_datafile(mar~ "${CMAKE_CURRENT_SOURCE_DIR}/mar~-help.pd")

file(GLOB audio_files "${CMAKE_CURRENT_SOURCE_DIR}/resources/*")
pd_add_datafile(mar~ "${audio_files}")
target_compile_options(mar_tilde PRIVATE -ggdb)
