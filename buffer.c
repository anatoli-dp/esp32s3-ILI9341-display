#include <esp_heap_caps.h>

#include "gc.h"

gc_color_t* display_buffer;
uint16_t* converted_buffer;

void buffer_init ()
{
    display_buffer = heap_caps_malloc(GC_DISPLAY_WIDTH * GC_DISPLAY_HEIGHT * sizeof(gc_color_t), MALLOC_CAP_SPIRAM);
    converted_buffer = heap_caps_malloc(GC_DISPLAY_HEIGHT * GC_DISPLAY_WIDTH * sizeof(uint16_t), MALLOC_CAP_SPIRAM);
}

static uint16_t convert_pixel (gc_color_t pixel)
{
    uint16_t r = ((pixel.r >> 3) & 0x1f) << 11;
    uint16_t g = ((pixel.g >> 2) & 0x3f) << 5;
    uint16_t b = (pixel.b >> 3) & 0x1f;

    return (uint16_t) (r | g | b);
}

void buffer_set_px (int posX, int posY, gc_color_t color)
{
    display_buffer[posY * GC_DISPLAY_WIDTH + posX] = color;
    converted_buffer[posY * GC_DISPLAY_WIDTH + posX] = convert_pixel(color);
}

gc_color_t buffer_get_px (int posX, int posY)
{
    return display_buffer[posY * GC_DISPLAY_WIDTH + posX];
}

void* buffer_get_converted ()
{
    /*for (int x = 0; x < GC_DISPLAY_WIDTH; x++)
    {
        for (int y = 0; y < GC_DISPLAY_HEIGHT; y++)
        {
            gc_color_t pixel = buffer_get_px(x, y);
            converted_buffer[y * GC_DISPLAY_WIDTH + x] = convert_pixel(pixel);
        }
    }*/

    return converted_buffer;
}