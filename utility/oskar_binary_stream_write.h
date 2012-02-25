/*
 * Copyright (c) 2011, The University of Oxford
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the University of Oxford nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef OSKAR_BINARY_STREAM_WRITE_H_
#define OSKAR_BINARY_STREAM_WRITE_H_

/**
 * @file oskar_binary_stream_write.h
 */

#include "oskar_global.h"

#ifdef __cplusplus
#include <cstdio>
#else
#include <stdio.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Writes a block of binary data to an output stream.
 *
 * @details
 * This function writes a block of binary data to an output stream.
 * The data are written at the current position of the stream, in native
 * byte order.
 *
 * The tag is specified as an extended tag, using a group name and a tag name
 * that are both given as strings.
 *
 * @param[in,out] stream   An output stream.
 * @param[in] data_type    Type of the memory (as in oskar_Mem).
 * @param[in] name_group   Tag group name.
 * @param[in] name_tag     Tag name.
 * @param[in] user_index   User-defined index.
 * @param[in] data_size    Size of data block, in bytes.
 * @param[out] data        Pointer to memory block to write.
 */
OSKAR_EXPORT
int oskar_binary_stream_write(FILE* stream, unsigned char data_type,
        const char* name_group, const char* name_tag, int user_index,
        size_t data_size, const void* data);

/**
 * @brief Writes a single double-precision value to an output stream.
 *
 * @details
 * This function writes a single double-precision value to an output stream.
 * The data are written at the current position of the stream, in native
 * byte order.
 *
 * The tag is specified as an extended tag, using a group name and a tag name
 * that are both given as strings.
 *
 * @param[in,out] stream   An output stream.
 * @param[in] name_group   Tag group name.
 * @param[in] name_tag     Tag name.
 * @param[in] user_index   User-defined index.
 * @param[in] value        Value to write.
 */
OSKAR_EXPORT
int oskar_binary_stream_write_double(FILE* stream, const char* name_group,
        const char* name_tag, int user_index, double value);

/**
 * @brief Writes a single integer value to an output stream.
 *
 * @details
 * This function writes a single integer value to an output stream.
 * The data are written at the current position of the stream, in native
 * byte order.
 *
 * The tag is specified as an extended tag, using a group name and a tag name
 * that are both given as strings.
 *
 * @param[in,out] stream   An output stream.
 * @param[in] name_group   Tag group name.
 * @param[in] name_tag     Tag name.
 * @param[in] user_index   User-defined index.
 * @param[in] value        Value to write.
 */
OSKAR_EXPORT
int oskar_binary_stream_write_int(FILE* stream, const char* name_group,
        const char* name_tag, int user_index, int value);

/**
 * @brief Writes a block of binary data to an output stream.
 *
 * @details
 * This function writes a block of binary data to an output stream.
 * The data are written at the current position of the stream, in native
 * byte order.
 *
 * The tag is specified as a standard tag, using a group ID and a tag ID
 * that are both given as bytes.
 *
 * @param[in,out] stream   An output stream.
 * @param[in] data_type    Type of the memory (as in oskar_Mem).
 * @param[in] id_group     Tag group identifier.
 * @param[in] id_tag       Tag identifier.
 * @param[in] user_index   User-defined index.
 * @param[in] data_size    Size of data block, in bytes.
 * @param[out] data        Pointer to memory block to write.
 */
OSKAR_EXPORT
int oskar_binary_stream_write_std(FILE* stream, unsigned char data_type,
        unsigned char id_group, unsigned char id_tag, int user_index,
        size_t data_size, const void* data);

/**
 * @brief Writes a single double-precision value to an output stream.
 *
 * @details
 * This function writes a single double-precision value to an output stream.
 * The data are written at the current position of the stream, in native
 * byte order.
 *
 * The tag is specified as a standard tag, using a group ID and a tag ID
 * that are both given as bytes.
 *
 * @param[in,out] stream   An output stream.
 * @param[in] id_group     Tag group identifier.
 * @param[in] id_tag       Tag identifier.
 * @param[in] user_index   User-defined index.
 * @param[in] value        Value to write.
 */
OSKAR_EXPORT
int oskar_binary_stream_write_std_double(FILE* stream, unsigned char id_group,
        unsigned char id_tag, int user_index, double value);

/**
 * @brief Writes a single integer value to an output stream.
 *
 * @details
 * This function writes a single integer value to an output stream.
 * The data are written at the current position of the stream, in native
 * byte order.
 *
 * The tag is specified as a standard tag, using a group ID and a tag ID
 * that are both given as bytes.
 *
 * @param[in,out] stream   An output stream.
 * @param[in] id_group     Tag group identifier.
 * @param[in] id_tag       Tag identifier.
 * @param[in] user_index   User-defined index.
 * @param[in] value        Value to write.
 */
OSKAR_EXPORT
int oskar_binary_stream_write_std_int(FILE* stream, unsigned char id_group,
        unsigned char id_tag, int user_index, int value);

#ifdef __cplusplus
}
#endif

#endif /* OSKAR_BINARY_STREAM_WRITE_H_ */