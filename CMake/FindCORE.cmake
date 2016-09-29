# Copyright (c) 2016, Blue Brain Project
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.


# FindCORE
# -------------
#
# Find CORE
#
# Find the CORE Blue Brain HPC utils library
#
# Using CORE:
#
# ::
#
#   find_package(CORE REQUIRED)
#   include_directories(${CORE_INCLUDE_DIRS})
#   target_link_libraries(foo ${CORE_LIBRARIES})
#
# This module sets the following variables:
#
# ::
#
#   CORE_FOUND - set to true if the library is found
#   CORE_INCLUDE_DIRS - list of required include directories
#   CORE_LIBRARIES - list of libraries to be linked

#=============================================================================
# Copyright 2015 Adrien Devresse <adrien.devresse@epfl.ch>
#
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)


# UNIX paths are standard, no need to write.
find_path(CORE_INCLUDE_DIR apf.h)

find_library(CORE_LIBRARY_APF NAMES apf)
find_library(CORE_LIBRARY_APF_ZOLTAN NAMES apf_zoltan)
find_library(CORE_LIBRARY_GMI NAMES gmi)
find_library(CORE_LIBRARY_LION NAMES lion)
find_library(CORE_LIBRARY_MA NAMES ma)
find_library(CORE_LIBRARY_MDS NAMES mds)
find_library(CORE_LIBRARY_MTH NAMES mth)
find_library(CORE_LIBRARY_PARM NAMES parma)
find_library(CORE_LIBRARY_PCU NAMES pcu)

#get_filename_component(CORE_LIB_DIR ${CORE_LIBRARY} DIRECTORY)

# Checks 'REQUIRED', 'QUIET' and versions.
include(FindPackageHandleStandardArgs)
list(APPEND CORE_LIBRARIES ${CORE_LIBRARY_PCU}
    ${CORE_LIBRARY_GMI} ${CORE_LIBRARY_MDS} ${CORE_LIBRARY_APF}
    ${CORE_LIBRARY_APF_ZOLTAN} ${CORE_LIBRARY_MA} ${CORE_LIBRARY_PARMA}
    ${CORE_LIBRARY_LION} ${CORE_LIBRARY_MTH})

find_package_handle_standard_args(CORE
  FOUND_VAR CORE_FOUND
  REQUIRED_VARS CORE_INCLUDE_DIR CORE_LIBRARIES)
  
