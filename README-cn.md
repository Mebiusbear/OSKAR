# pierce_points

## oskar_evaluate_tec_tid
```C++
/* TID */
struct oskar_SettingsTIDscreen
{
    double height_km;       /* Height of the TID screen, km */
    int num_components;     /* Number of TID components in the screen */
    double* amp;            /* Relative amplitude compared to TEC0 */
    double* wavelength;     /* km */
    double* speed;          /* km/h */
    double* theta;          /* deg. */
};
```
```C++
amp = TID->amp[i];
w = TID->wavelength[i] / (earth_radius + TID->height_km);
th = TID->theta[i] * M_PI/180.;
v = (TID->speed[i]/(earth_radius + TID->height_km)) / 3600;

pp_lon = oskar_mem_double_const(lon, &status)[j];
pp_lat = oskar_mem_double_const(lat, &status)[j];
pp_sec = oskar_mem_double_const(rel_path_length, &status)[j];
```

```C++
struct oskar_WorkJonesZ
{
    oskar_Mem* hor_x;        /* Pierce point horizontal x coordinate */
    oskar_Mem* hor_y;        /* Pierce point horizontal y coordinate */
    oskar_Mem* hor_z;        /* Pierce point horizontal z coordinate */

    oskar_Mem* pp_lon;       /* Pierce point longitude, in radians */
    oskar_Mem* pp_lat;       /* Pierce point latitude, in radians */
    oskar_Mem* pp_rel_path;  /* Pierce point relative path length.
                                (the extra path, relative to the vertical for
                                the ionospheric column defined by the pierce
                                point) */

    oskar_Mem* screen_TEC;    /* TEC screen values for each pierce point */
    oskar_Mem* total_TEC;     /* Total TEC values for each pierce point */
};

/* Evaluate TEC values for the screen */
oskar_evaluate_tec_tid(work->screen_TEC, num_pp, work->pp_lon, work->pp_lat, work->pp_rel_path, settings->TEC0, &settings->TID[i], gast);

oskar_mem_add(work->total_TEC, work->total_TEC, work->screen_TEC, 0, 0, 0, oskar_mem_length(work->total_TEC), status);
```


---





```c++
/* function(void oskar_station_beam) */
tec_phase = oskar_station_work_evaluate_tec_screen(work,
        (int) num_points_orig, enu[0], enu[1], u, v,
        time_index, frequency_hz, status);
```

```c++
/* function(const oskar_Mem* oskar_station_work_evaluate_tec_screen) */
oskar_evaluate_tec_screen(work->isoplanatic_screen,
        num_points, l, m, station_u_m, station_v_m,
        frequency_hz, work->screen_height_km * 1000.0,
        work->screen_pixel_size_m,
        work->screen_num_pixels_x, work->screen_num_pixels_y,
        work->tec_screen, 0, work->screen_output, status);
```

```c++
/* function(void oskar_evaluate_tec_screen) */
evaluate_tec_screen_double(isoplanatic, num_points,
        oskar_mem_double_const(l, status),
        oskar_mem_double_const(m, status),
        station_u_m, station_v_m, inv_freq_hz,
        screen_height_m, inv_pixel_size_m,
        screen_num_pixels_x, screen_num_pixels_y,
        oskar_mem_double_const(tec_screen, status), offset_out,
        oskar_mem_double2(out, status));
```

```c++
/* KERNEL */
#define OSKAR_EVALUATE_TEC_SCREEN(NAME, FP, FP2)\
KERNEL(NAME) (const int isoplanatic,\
        const int num_points, GLOBAL_IN(FP, l), GLOBAL_IN(FP, m),\
        const FP station_u, const FP station_v, const FP inv_frequency_hz,\
        const FP screen_height_m, const FP inv_pixel_size_m,\
        const int screen_num_pixels_x, const int screen_num_pixels_y,\
        GLOBAL_IN(FP, screen), const int offset_out, GLOBAL_OUT(FP2, out))\
{\
    const int screen_half_x = screen_num_pixels_x / 2;\
    const int screen_half_y = screen_num_pixels_y / 2;\
    KERNEL_LOOP_X(int, i, 0, num_points)\
    FP2 comp;\
    FP s_l, s_m;\
    if (isoplanatic) {\
        s_l = (FP)0;\
        s_m = (FP)0;\
    }\
    else {\
        s_l = l[i];\
        s_m = m[i];\
    }\
    const FP world_x = (station_u + s_l * screen_height_m) * inv_pixel_size_m;\
    const FP world_y = (station_v + s_m * screen_height_m) * inv_pixel_size_m;\
    const int pix_x = screen_half_x + ROUND(FP, world_x);\
    const int pix_y = screen_half_y + ROUND(FP, world_y);\
    if (pix_x >= 0 && pix_y >= 0 &&\
            pix_x < screen_num_pixels_x && pix_y < screen_num_pixels_y)\
    {\
        FP sin_phase, cos_phase;\
        const FP tec = screen[pix_x + pix_y * screen_num_pixels_x];\
        const FP phase = tec * ((FP) -8.44797245e9) * inv_frequency_hz;\
        SINCOS(phase, sin_phase, cos_phase);\
        comp.x = cos_phase; comp.y = sin_phase;\
    }\
    else\
    {\
        comp.x = (FP)1; comp.y = (FP)0;\
    }\
    out[i + offset_out] = comp;\
    KERNEL_LOOP_END\
}\
OSKAR_REGISTER_KERNEL(NAME)
```

```c++
void oskar_mem_read_fits(oskar_Mem* data, size_t offset, size_t num_pixels,
        const char* file_name, int num_index_dims, const int* start_index,
        int* num_axes, int** axis_size, double** axis_inc, int* status)
```