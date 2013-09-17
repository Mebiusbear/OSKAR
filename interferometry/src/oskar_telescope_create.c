/*
 * Copyright (c) 2013, The University of Oxford
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

#include <private_telescope.h>
#include <oskar_telescope.h>

#include <oskar_mem.h>

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

oskar_Telescope* oskar_telescope_create(int type, int location,
        int num_stations, int* status)
{
    int i = 0;
    oskar_Telescope* telescope;

    /* Check all inputs. */
    if (!status)
    {
        oskar_set_invalid_argument(status);
        return 0;
    }

    /* Check type and location. */
    if (type != OSKAR_SINGLE && type != OSKAR_DOUBLE)
    {
        *status = OSKAR_ERR_BAD_DATA_TYPE;
        return 0;
    }
    if (location != OSKAR_LOCATION_CPU && location != OSKAR_LOCATION_GPU)
    {
        *status = OSKAR_ERR_BAD_LOCATION;
        return 0;
    }

    /* Allocate the data structure. */
    telescope = (oskar_Telescope*) malloc(sizeof(oskar_Telescope));
    if (!telescope)
    {
        *status = OSKAR_ERR_MEMORY_ALLOC_FAILURE;
        return 0;
    }

    /* Initialise private meta-data. */
    telescope->precision = type;
    telescope->location = location;

    /* Initialise the meta-data. */
    telescope->num_stations = num_stations;
    telescope->max_station_size = 0;
    telescope->max_station_depth = 1;
    telescope->coord_units = OSKAR_METRES;
    telescope->identical_stations = 0;
    telescope->use_common_sky = 1;
    telescope->seed_time_variable_station_element_errors = 0;
    telescope->longitude_rad = 0.0;
    telescope->latitude_rad = 0.0;
    telescope->altitude_m = 0.0;
    telescope->ra0_rad = 0.0;
    telescope->dec0_rad = 0.0;
    telescope->wavelength_m = 0.0;
    telescope->bandwidth_hz = 0.0;
    telescope->time_average_sec = 0.0;

    /* Initialise the arrays. */
    oskar_mem_init(&telescope->station_x, type, location,
            num_stations, 1, status);
    oskar_mem_init(&telescope->station_y, type, location,
            num_stations, 1, status);
    oskar_mem_init(&telescope->station_z, type, location,
            num_stations, 1, status);
    oskar_mem_init(&telescope->station_x_hor, type, location,
            num_stations, 1, status);
    oskar_mem_init(&telescope->station_y_hor, type, location,
            num_stations, 1, status);
    oskar_mem_init(&telescope->station_z_hor, type, location,
            num_stations, 1, status);

    /* Initialise the station structures. */
    telescope->station = NULL;
    if (num_stations > 0)
        telescope->station = malloc(num_stations * sizeof(oskar_Station*));
    for (i = 0; i < num_stations; ++i)
    {
        telescope->station[i] = oskar_station_create(type, location, 0, status);
    }

    /* Return pointer to data structure. */
    return telescope;
}

#ifdef __cplusplus
}
#endif