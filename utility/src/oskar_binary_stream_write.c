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

#include "utility/oskar_BinaryTag.h"
#include "utility/oskar_binary_stream_write.h"
#include "utility/oskar_endian.h"
#include "utility/oskar_Mem.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

int oskar_binary_stream_write(FILE* stream, unsigned char data_type,
        const char* name_group, const char* name_tag, int user_index,
        size_t data_size, const void* data)
{
    oskar_BinaryTag tag;
    size_t block_size;
    int tag_index;

    /* Initialise the tag. */
    char magic[] = "TAG";
    strcpy(tag.magic, magic);
    memset(tag.size_bytes, 0, sizeof(tag.size_bytes));
    memset(tag.user_index, 0, sizeof(tag.user_index));

    /* Set up the tag identifiers */
    tag.flags = (unsigned char)1;
    tag.data_type = data_type;
    tag.group.bytes = 1 + strlen(name_group);
    tag.tag.bytes = 1 + strlen(name_tag);

    /* Get the number of bytes in the block in little-endian byte order. */
    block_size = data_size + tag.group.bytes + tag.tag.bytes;
    if (sizeof(size_t) != 4 && sizeof(size_t) != 8)
    {
        return OSKAR_ERR_BAD_BINARY_FORMAT;
    }
    if (oskar_endian() != OSKAR_LITTLE_ENDIAN)
    {
        if (sizeof(size_t) == 4)
            oskar_endian_swap_4((char*)&block_size);
        else if (sizeof(size_t) == 8)
            oskar_endian_swap_8((char*)&block_size);
    }

    /* Copy the block size in bytes to the tag, as little endian. */
    memcpy(tag.size_bytes, &block_size, sizeof(size_t));

    /* Get the user index in little-endian byte order. */
    tag_index = user_index;
    if (oskar_endian() != OSKAR_LITTLE_ENDIAN)
    {
        if (sizeof(int) == 2)
            oskar_endian_swap_2((char*)&tag_index);
        else if (sizeof(int) == 4)
            oskar_endian_swap_4((char*)&tag_index);
        else if (sizeof(int) == 8)
            oskar_endian_swap_8((char*)&tag_index);
    }

    /* Copy the user index to the tag, as little endian. */
    memcpy(tag.user_index, &tag_index, sizeof(int));

    /* Write the tag to the file. */
    if (fwrite(&tag, sizeof(oskar_BinaryTag), 1, stream) != 1)
        return OSKAR_ERR_FILE_IO;

    /* Write the group name and tag name to the file. */
    if (fwrite(name_group, 1, tag.group.bytes, stream) != tag.group.bytes)
        return OSKAR_ERR_FILE_IO;
    if (fwrite(name_tag, 1, tag.tag.bytes, stream) != tag.tag.bytes)
        return OSKAR_ERR_FILE_IO;

    /* Write the data to the file. */
    if (fwrite(data, 1, data_size, stream) != data_size)
        return OSKAR_ERR_FILE_IO;

    return OSKAR_SUCCESS;
}

int oskar_binary_stream_write_double(FILE* stream, const char* name_group,
        const char* name_tag, int user_index, double value)
{
    return oskar_binary_stream_write(stream, OSKAR_DOUBLE, name_group,
            name_tag, user_index, sizeof(double), &value);
}

int oskar_binary_stream_write_int(FILE* stream, const char* name_group,
        const char* name_tag, int user_index, int value)
{
    return oskar_binary_stream_write(stream, OSKAR_INT, name_group,
            name_tag, user_index, sizeof(int), &value);
}

int oskar_binary_stream_write_std(FILE* stream, unsigned char data_type,
        unsigned char id_group, unsigned char id_tag, int user_index,
        size_t data_size, const void* data)
{
    oskar_BinaryTag tag;
    size_t block_size;
    int tag_index;

    /* Initialise the tag. */
    char magic[] = "TAG";
    strcpy(tag.magic, magic);
    memset(tag.size_bytes, 0, sizeof(tag.size_bytes));
    memset(tag.user_index, 0, sizeof(tag.user_index));

    /* Set up the tag identifiers */
    tag.flags = 0;
    tag.data_type = data_type;
    tag.group.id = id_group;
    tag.tag.id = id_tag;

    /* Get the number of bytes in the block in little-endian byte order. */
    block_size = data_size;
    if (sizeof(size_t) != 4 && sizeof(size_t) != 8)
    {
        return OSKAR_ERR_BAD_BINARY_FORMAT;
    }
    if (oskar_endian() != OSKAR_LITTLE_ENDIAN)
    {
        if (sizeof(size_t) == 4)
            oskar_endian_swap_4((char*)&block_size);
        else if (sizeof(size_t) == 8)
            oskar_endian_swap_8((char*)&block_size);
    }

    /* Copy the block size in bytes to the tag, as little endian. */
    memcpy(tag.size_bytes, &block_size, sizeof(size_t));

    /* Get the user index in little-endian byte order. */
    tag_index = user_index;
    if (oskar_endian() != OSKAR_LITTLE_ENDIAN)
    {
        if (sizeof(int) == 2)
            oskar_endian_swap_2((char*)&tag_index);
        else if (sizeof(int) == 4)
            oskar_endian_swap_4((char*)&tag_index);
        else if (sizeof(int) == 8)
            oskar_endian_swap_8((char*)&tag_index);
    }

    /* Copy the user index to the tag, as little endian. */
    memcpy(tag.user_index, &tag_index, sizeof(int));

    /* Write the tag to the file. */
    if (fwrite(&tag, sizeof(oskar_BinaryTag), 1, stream) != 1)
        return OSKAR_ERR_FILE_IO;

    /* Write the data to the file. */
    if (fwrite(data, 1, data_size, stream) != data_size)
        return OSKAR_ERR_FILE_IO;

    return OSKAR_SUCCESS;
}

int oskar_binary_stream_write_std_double(FILE* stream, unsigned char id_group,
        unsigned char id_tag, int user_index, double value)
{
    return oskar_binary_stream_write_std(stream, OSKAR_DOUBLE, id_group, id_tag,
            user_index, sizeof(double), &value);
}

int oskar_binary_stream_write_std_int(FILE* stream, unsigned char id_group,
        unsigned char id_tag, int user_index, int value)
{
    return oskar_binary_stream_write_std(stream, OSKAR_INT, id_group, id_tag,
            user_index, sizeof(int), &value);
}

#ifdef __cplusplus
}
#endif