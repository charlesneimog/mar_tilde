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
    file(DOWNLOAD "https://github.com/pure-data/pd.cmake/releases/download/v2.0.11/pd.cmake" ${PDCMAKE_FILE})
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
cpmaddpackage("gh:jprjr/miniflac@1.1.1")
if(miniflac_ADDED)
    include_directories("${miniflac_SOURCE_DIR}")
endif()

# ╭──────────────────────────────────────╮
# │              Resampler               │
# ╰──────────────────────────────────────╯
cpmaddpackage("gh:avaneev/r8brain-free-src#master")
if(r8brain-free-src_ADDED)
    include_directories("${r8brain-free-src_SOURCE_DIR}")
endif()

# ╭──────────────────────────────────────╮
# │              Pd Object               │
# ╰──────────────────────────────────────╯
pd_add_external(mar~ "${CMAKE_CURRENT_SOURCE_DIR}/mar~.cpp" LINK_LIBRARIES AudioFile)
pd_add_datafile(mar~ "${CMAKE_CURRENT_SOURCE_DIR}/mar~-help.pd")

file(GLOB audio_files "${CMAKE_CURRENT_SOURCE_DIR}/resources/*")
pd_add_datafile(mar~ "${audio_files}")
