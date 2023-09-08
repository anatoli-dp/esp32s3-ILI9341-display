#include <esp_log.h>

#include <driver/ledc.h>
#include <driver/spi_master.h>
#include <esp_lcd_panel_io.h>
#include <esp_lcd_panel_ops.h>

#include "gc.h"

#include "ILI9341.h"
#include "buffer.h"
#include "display.h"

static const char* TAG = "DISPLAY";

#define GPIO_MOSI 10
#define GPIO_CLOCK 11
#define GPIO_CS 3
#define GPIO_RESET 46
#define GPIO_DC 9
#define GPIO_BACKLIGHT 12

static const int DISPLAY_SPI_QUEUE_LEN = 10;
static const int SPI_MAX_TRANSFER_SIZE = 32768;

static const ledc_mode_t BACKLIGHT_LEDC_MODE = LEDC_LOW_SPEED_MODE;
static const ledc_channel_t BACKLIGHT_LEDC_CHANNEL = LEDC_CHANNEL_0;
static const ledc_timer_t BACKLIGHT_LEDC_TIMER = LEDC_TIMER_1;
static const ledc_timer_bit_t BACKLIGHT_LEDC_TIMER_RESOLUTION = LEDC_TIMER_10_BIT;
static const uint32_t BACKLIGHT_LEDC_FRQUENCY = 5000;

static esp_lcd_panel_io_handle_t lcd_io_handle = NULL;
static esp_lcd_panel_handle_t lcd_handle = NULL;

void display_init ()
{
    //initialize brightness with ledc
    const ledc_channel_config_t LCD_backlight_channel = {
            .gpio_num = (gpio_num_t)GPIO_BACKLIGHT,
            .speed_mode = BACKLIGHT_LEDC_MODE,
            .channel = BACKLIGHT_LEDC_CHANNEL,
            .intr_type = LEDC_INTR_DISABLE,
            .timer_sel = BACKLIGHT_LEDC_TIMER,
            .duty = 0,
            .hpoint = 0,
            .flags = {
                    .output_invert = 0
            }
    };
    const ledc_timer_config_t LCD_backlight_timer = {
            .speed_mode = BACKLIGHT_LEDC_MODE,
            .duty_resolution = BACKLIGHT_LEDC_TIMER_RESOLUTION,
            .timer_num = BACKLIGHT_LEDC_TIMER,
            .freq_hz = BACKLIGHT_LEDC_FRQUENCY,
            .clk_cfg = LEDC_AUTO_CLK
    };
    ESP_LOGI(TAG, "Initializing LEDC for backlight pin: %d", GPIO_BACKLIGHT);

    ESP_ERROR_CHECK(ledc_timer_config(&LCD_backlight_timer));
    ESP_ERROR_CHECK(ledc_channel_config(&LCD_backlight_channel));

    gc_display_set_brightness(0);

    //initialize spi
    ESP_LOGI(TAG, "Initializing SPI bus (MOSI:%d, CLK:%d)", GPIO_MOSI, GPIO_CLOCK);
    spi_bus_config_t bus = {
            .mosi_io_num = GPIO_MOSI,
            .sclk_io_num = GPIO_CLOCK,
            .quadwp_io_num = GPIO_NUM_NC,
            .quadhd_io_num = GPIO_NUM_NC,
            .data4_io_num = GPIO_NUM_NC,
            .data5_io_num = GPIO_NUM_NC,
            .data6_io_num = GPIO_NUM_NC,
            .data7_io_num = GPIO_NUM_NC,
            .max_transfer_sz = SPI_MAX_TRANSFER_SIZE,
            .flags = SPICOMMON_BUSFLAG_SCLK | SPICOMMON_BUSFLAG_MISO | SPICOMMON_BUSFLAG_MOSI | SPICOMMON_BUSFLAG_MASTER,
            .intr_flags = ESP_INTR_FLAG_LOWMED | ESP_INTR_FLAG_IRAM
    };

    ESP_ERROR_CHECK(spi_bus_initialize(SPI3_HOST, &bus, SPI_DMA_CH_AUTO));

    //initialize display
    const esp_lcd_panel_io_spi_config_t io_config = {
            .cs_gpio_num = GPIO_CS,
            .dc_gpio_num = GPIO_DC,
            .spi_mode = 0,
            .pclk_hz = 60000000,
            .trans_queue_depth = DISPLAY_SPI_QUEUE_LEN,
            .lcd_cmd_bits = 8,
            .lcd_param_bits = 8,
            .flags = {
                    .dc_low_on_data = 0,
                    .octal_mode = 0,
                    .sio_mode = 0,
                    .lsb_first = 0,
                    .cs_high_active = 0
            }
    };

    const esp_lcd_panel_dev_config_t lcd_config = {
            .reset_gpio_num = GPIO_RESET,
            .color_space = LCD_COLOR_SPACE_RGB,
            .bits_per_pixel = 16,
            .flags = {
                    .reset_active_high = 0
            },
            .vendor_config = NULL
    };

    ESP_ERROR_CHECK(esp_lcd_new_panel_io_spi((esp_lcd_spi_bus_handle_t)SPI3_HOST, &io_config, &lcd_io_handle));

    ESP_ERROR_CHECK(esp_lcd_new_panel_ili9341(lcd_io_handle, &lcd_config, &lcd_handle));

    ESP_ERROR_CHECK(esp_lcd_panel_reset(lcd_handle));
    ESP_ERROR_CHECK(esp_lcd_panel_init(lcd_handle));
    ESP_ERROR_CHECK(esp_lcd_panel_invert_color(lcd_handle, false));
    ESP_ERROR_CHECK(esp_lcd_panel_swap_xy(lcd_handle, true));
    ESP_ERROR_CHECK(esp_lcd_panel_mirror(lcd_handle, false, true));
    ESP_ERROR_CHECK(esp_lcd_panel_set_gap(lcd_handle, 0, 0));
    ESP_ERROR_CHECK(esp_lcd_panel_disp_on_off(lcd_handle, true));

    buffer_init();

    gc_display_set_brightness(100);
}

void gc_display_update ()
{
    uint16_t* buffer = buffer_get_converted();
    esp_lcd_panel_draw_bitmap(lcd_handle, 0, 0, 320, 240, buffer);
}

void gc_display_set_brightness (int brightness_percentage)
{
    if (brightness_percentage > 100)
    {
        brightness_percentage = 100;
    }
    else if (brightness_percentage < 0)
    {
        brightness_percentage = 0;
    }
    ESP_LOGI(TAG, "Setting backlight to %d%%", brightness_percentage);

    // LEDC resolution set to 10bits, thus: 100% = 1023
    uint32_t duty_cycle = (1023 * brightness_percentage) / 100;
    ESP_ERROR_CHECK(ledc_set_duty(BACKLIGHT_LEDC_MODE, BACKLIGHT_LEDC_CHANNEL, duty_cycle));
    ESP_ERROR_CHECK(ledc_update_duty(BACKLIGHT_LEDC_MODE, BACKLIGHT_LEDC_CHANNEL));
}