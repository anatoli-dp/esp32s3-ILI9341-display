#pragma once

#include "gc.h"

void buffer_init();

void buffer_set_px(int posX, int posY, gc_color_t color);
gc_color_t buffer_get_px(int posX, int posY);

void* buffer_get_converted();