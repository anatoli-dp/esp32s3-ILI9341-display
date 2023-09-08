#include <math.h>

#include <esp_heap_caps.h>

#include "gc.h"
#include "buffer.h"

int abs (int n) {return ((n > 0) ? n : (n * (-1)));}

gc_color_t currentColor;

void gc_set_color (gc_color_t color)
{
    currentColor = color;
}

void gc_draw_px (int posX, int posY)
{
    if (posX >= 0 && posX <= GC_DISPLAY_WIDTH)
    {
        if (posY >= 0 && posY <= GC_DISPLAY_HEIGHT)
        {
            buffer_set_px(posX, posY, currentColor);
        }
    }
}

void gc_draw_line (int x0, int y0, int x1, int y1)
{
    int dy = y1 - y0;
    int dx = x1 - x0;
    int stepx, stepy;

    if (dy < 0) { dy = -dy; stepy = -1; } else { stepy = 1; }
    if (dx < 0) { dx = -dx; stepx = -1; } else { stepx = 1; }

    dy <<= 1;
    dx <<= 1;
    gc_draw_px(x0, y0);

    if (dx > dy) {
        int fraction = dy - (dx >> 1);
        while (x0 != x1)
        {
            if (fraction >= 0)
            {
                y0 += stepy;
                fraction -= dx;
            }
            x0 += stepx;
            fraction += dy;
            gc_draw_px(x0, y0);
        }
    } else {
        int fraction = dx - (dy >> 1);
        while (y0 != y1)
        {
            if (fraction >= 0)
            {
                x0 += stepx;
                fraction -= dy;
            }
            y0 += stepy;
            fraction += dx;
            gc_draw_px(x0, y0);
        }
    }
}

void gc_draw_rect (int x0, int y0, int x1, int y1)
{
    gc_draw_line(x0, y0, x1, y0);
    gc_draw_line(x1, y0, x1, y1);
    gc_draw_line(x1, y1, x0, y1);
    gc_draw_line(x0, y1, x0, y0);
}

void gc_draw_rect_filled (int x0, int y0, int x1, int y1)
{
    gc_draw_rect(x0, y0, x1, y1);
    int x = round((x0 + x1) / 2);
    int y = round((y0 + y1) / 2);
    gc_draw_fill(x, y);
}

void gc_draw_circle (int xm, int ym, int r)
{
    int x = -r, y = 0, err = 2-2*r;
    do {
        gc_draw_px(xm-x, ym+y);
        gc_draw_px(xm-y, ym-x);
        gc_draw_px(xm+x, ym-y);
        gc_draw_px(xm+y, ym+x);
        r = err;
        if (r <= y) err += ++y*2+1;
        if (r > x || err > y)
            err += ++x*2+1;
    } while (x < 0);
}

void gc_draw_circle_filled (int xm, int ym, int r)
{
    gc_draw_circle(xm, ym, r);
    gc_draw_fill(xm, ym);
}

void gc_draw_ellipse (int xm, int ym, int a, int b)
{
    long x = -a, y = 0;
    long e2 = b, dx = (1+2*x)*e2*e2;
    long dy = x*x, err = dx+dy;

    do {
        gc_draw_px(xm-x, ym+y);
        gc_draw_px(xm+x, ym+y);
        gc_draw_px(xm+x, ym-y);
        gc_draw_px(xm-x, ym-y);
        e2 = 2*err;
        if (e2 >= dx) { x++; err += dx += 2*(long)b*b; }
        if (e2 <= dy) { y++; err += dy += 2*(long)a*a; }
    } while (x <= 0);

    while (y++ < b)
    {
        gc_draw_px(xm, ym+y);
        gc_draw_px(xm, ym-y);
    }
}

void gc_draw_ellipse_filled (int xm, int ym, int a, int b)
{
    gc_draw_ellipse(xm, ym, a, b);
    gc_draw_fill(xm, ym);
}

static void plotEllipseRect(int x0, int y0, int x1, int y1)
{
    long a = abs(x1-x0), b = abs(y1-y0), b1 = b&1;
    double dx = 4*(1.0-a)*b*b, dy = 4*(b1+1)*a*a;
    double err = dx+dy+b1*a*a, e2;

    if (x0 > x1) { x0 = x1; x1 += a; }
    if (y0 > y1) y0 = y1;
    y0 += (b+1)/2; y1 = y0-b1;
    a = 8*a*a; b1 = 8*b*b;

    do {
        gc_draw_px(x1, y0);
        gc_draw_px(x0, y0);
        gc_draw_px(x0, y1);
        gc_draw_px(x1, y1);
        e2 = 2*err;
        if (e2 <= dy) { y0++; y1--; err += dy += a; }
        if (e2 >= dx || 2*err > dy) { x0++; x1--; err += dx += b1; }
    } while (x0 <= x1);

    while (y0-y1 <= b)
    {
        gc_draw_px(x0-1, y0);
        gc_draw_px(x1+1, y0++);
        gc_draw_px(x0-1, y1);
        gc_draw_px(x1+1, y1--);
    }
}

static void plotQuadRationalBezierSeg (int x0, int y0, int x1, int y1, int x2, int y2, float w)
{
    int sx = x2-x1, sy = y2-y1;
    double dx = x0-x2, dy = y0-y2, xx = x0-x1, yy = y0-y1;
    double xy = xx*sy+yy*sx, cur = xx*sy-yy*sx, err;

    assert(xx*sx <= 0.0 && yy*sy <= 0.0);

    if (cur != 0.0 && w > 0.0)
    {
        if (sx*(long)sx+sy*(long)sy > xx*xx+yy*yy)
        {
            x2 = x0; x0 -= dx; y2 = y0; y0 -= dy; cur = -cur;
        }
        
        xx = 2.0*(4.0*w*sx*xx+dx*dx);
        yy = 2.0*(4.0*w*sy*yy+dy*dy);
        sx = x0 < x2 ? 1 : -1;
        sy = y0 < y2 ? 1 : -1;
        xy = -2.0*sx*sy*(2.0*w*xy+dx*dy);

        if (cur*sx*sy < 0.0)
        {
            xx = -xx; yy = -yy; xy = -xy; cur = -cur;
        }
        
        dx = 4.0*w*(x1-x0)*sy*cur+xx/2.0+xy;
        dy = 4.0*w*(y0-y1)*sx*cur+yy/2.0+xy;

        if (w < 0.5 && (dy > xy || dx < xy))
        {
            cur = (w+1.0)/2.0; w = sqrt(w); xy = 1.0/(w+1.0);
            sx = floor((x0+2.0*w*x1+x2)*xy/2.0+0.5);
            sy = floor((y0+2.0*w*y1+y2)*xy/2.0+0.5);
            dx = floor((w*x1+x0)*xy+0.5); dy = floor((y1*w+y0)*xy+0.5);
            plotQuadRationalBezierSeg(x0,y0, dx,dy, sx,sy, cur);
            dx = floor((w*x1+x2)*xy+0.5); dy = floor((y1*w+y2)*xy+0.5);
            plotQuadRationalBezierSeg(sx,sy, dx,dy, x2,y2, cur);
            return;
        }
        err = dx+dy-xy;
        do {
            gc_draw_px(x0,y0);
            if (x0 == x2 && y0 == y2) return;
            x1 = 2*err > dy; y1 = 2*(err+yy) < -dy;
            if (2*err < dx || y1) { y0 += sy; dy += xy; err += dx += xx; }
            if (2*err > dx || x1) { x0 += sx; dx += xy; err += dy += yy; }
        } while (dy <= xy && dx >= xy);
    } 
    
    gc_draw_line(x0,y0, x2,y2);
}

static void plotRotatedEllipseRect(int x0, int y0, int x1, int y1, long zd)
{
    int xd = x1-x0, yd = y1-y0;
    float w = xd*(long)yd;
    if (zd == 0) return plotEllipseRect(x0,y0, x1,y1);
    if (w != 0.0) w = (w-zd)/(w+w);
    assert(w <= 1.0 && w >= 0.0);
    xd = floor(xd*w+0.5); yd = floor(yd*w+0.5);
    plotQuadRationalBezierSeg(x0,y0+yd, x0,y0, x0+xd,y0, 1.0-w);
    plotQuadRationalBezierSeg(x0,y0+yd, x0,y1, x1-xd,y1, w);
    plotQuadRationalBezierSeg(x1,y1-yd, x1,y1, x1-xd,y1, 1.0-w);
    plotQuadRationalBezierSeg(x1,y1-yd, x1,y0, x0+xd,y0, w);
}

void gc_draw_ellipse_rot (int x, int y, int a, int b, float angle)
{
    float xd = (long)a*a, yd = (long)b*b;
    float s = sin(angle), zd = (xd-yd)*s;
    xd = sqrt(xd-zd*s), yd = sqrt(yd+zd*s);
    a = xd+0.5; b = yd+0.5; zd = zd*a*b/(xd*yd);
    plotRotatedEllipseRect(x-a,y-b, x+a,y+b, (long)(4*zd*cos(angle)));
}

void gc_draw_ellipse_rot_filled (int x, int y, int a, int b, float angle)
{
    gc_draw_ellipse_rot(x, y, a, b, angle);
    gc_draw_fill(x, y);
}

static void quadBezierSeg (int x0, int y0, int x1, int y1, int x2, int y2)
{
    int sx = x2-x1, sy = y2-y1;
    long xx = x0-x1, yy = y0-y1, xy;
    double dx, dy, err, cur = xx*sy-yy*sx;

    assert(xx*sx <= 0 && yy*sy <= 0);

    if (sx*(long)sx+sy*(long)sy > xx*xx+yy*yy)
    {
        x2 = x0; x0 = sx+x1; y2 = y0; y0 = sy+y1; cur = -cur;
    }

    if (cur != 0)
    {
        xx += sx; xx *= sx = x0 < x2 ? 1 : -1;
        yy += sy; yy *= sy = y0 < y2 ? 1 : -1;
        xy = 2*xx*yy; xx *= xx; yy *= yy;

        if (cur*sx*sy < 0)
        {
            xx = -xx; yy = -yy; xy = -xy; cur = -cur;
        }

        dx = 4.0*sy*cur*(x1-x0)+xx-xy;
        dy = 4.0*sx*cur*(y0-y1)+yy-xy;
        xx += xx; yy += yy; err = dx+dy+xy;

        do {
            gc_draw_px(x0, y0);
            if (x0 == x2 && y0 == y2) return;
            y1 = 2*err < dx;
            if (2*err > dy) { x0 += sx; dx -= xy; err += dy += yy; }
            if (    y1    ) { y0 += sy; dy -= xy; err += dx += xx; }
        } while (dy < 0 && dx > 0);
    }

    gc_draw_line(x0, y0, x2, y2);
}

void gc_draw_quadratic_bezier (int x0, int y0, int x1, int y1, int x2, int y2)
{
    int x = x0-x1, y = y0-y1;
    double t = x0-2*x1+x2, r;

    if ((long)x*(x2-x1) > 0)
    {
        if ((long)y*(y2-y1) > 0)
            if (fabs((y0-2*y1+y2)/t*x) > abs(y)) {
                x0 = x2; x2 = x+x1; y0 = y2; y2 = y+y1;
            }
        t = (x0-x1)/t;
        r = (1-t)*((1-t)*y0+2.0*t*y1)+t*t*y2;
        t = (x0*x2-x1*x1)*t/(x0-x1);
        x = floor(t+0.5); y = floor(r+0.5);
        r = (y1-y0)*(t-x0)/(x1-x0)+y0;
        quadBezierSeg(x0,y0, x,floor(r+0.5), x,y);
        r = (y1-y2)*(t-x2)/(x1-x2)+y2;
        x0 = x1 = x; y0 = y; y1 = floor(r+0.5);
    }

    if ((long)(y0-y1)*(y2-y1) > 0)
    {
        t = y0-2*y1+y2; t = (y0-y1)/t;
        r = (1-t)*((1-t)*x0+2.0*t*x1)+t*t*x2;
        t = (y0*y2-y1*y1)*t/(y0-y1);
        x = floor(r+0.5); y = floor(t+0.5);
        r = (x1-x0)*(t-y0)/(y1-y0)+x0;
        quadBezierSeg(x0,y0, floor(r+0.5),y, x,y);
        r = (x1-x2)*(t-y2)/(y1-y2)+x2;
        x0 = x; x1 = floor(r+0.5); y0 = y1 = y;
    }

    quadBezierSeg(x0,y0, x1,y1, x2,y2);
}

static void cubicBezierSeg (int x0, int y0, float x1, float y1, float x2, float y2, int x3, int y3)
{
    int f, fx, fy, leg = 1;
    int sx = x0 < x3 ? 1 : -1, sy = y0 < y3 ? 1 : -1;
    float xc = -fabs(x0+x1-x2-x3), xa = xc-4*sx*(x1-x2), xb = sx*(x0-x1-x2+x3);
    float yc = -fabs(y0+y1-y2-y3), ya = yc-4*sy*(y1-y2), yb = sy*(y0-y1-y2+y3);
    double ab, ac, bc, cb, xx, xy, yy, dx, dy, ex, *pxy, EP = 0.01;

    assert((x1-x0)*(x2-x3) < EP && ((x3-x0)*(x1-x2) < EP || xb*xb < xa*xc+EP));
    assert((y1-y0)*(y2-y3) < EP && ((y3-y0)*(y1-y2) < EP || yb*yb < ya*yc+EP));

    if (xa == 0 && ya == 0)
     {
        sx = floor((3*x1-x0+1)/2); sy = floor((3*y1-y0+1)/2);
        return quadBezierSeg(x0,y0, sx,sy, x3,y3);
    }

    x1 = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+1;
    x2 = (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)+1;

    do {
        ab = xa*yb-xb*ya; ac = xa*yc-xc*ya; bc = xb*yc-xc*yb;
        ex = ab*(ab+ac-3*bc)+ac*ac;
        f = ex > 0 ? 1 : sqrt(1+1024/x1);
        ab *= f; ac *= f; bc *= f; ex *= f*f;
        xy = 9*(ab+ac+bc)/8; cb = 8*(xa-ya);
        dx = 27*(8*ab*(yb*yb-ya*yc)+ex*(ya+2*yb+yc))/64-ya*ya*(xy-ya);
        dy = 27*(8*ab*(xb*xb-xa*xc)-ex*(xa+2*xb+xc))/64-xa*xa*(xy+xa);
        xx = 3*(3*ab*(3*yb*yb-ya*ya-2*ya*yc)-ya*(3*ac*(ya+yb)+ya*cb))/4;
        yy = 3*(3*ab*(3*xb*xb-xa*xa-2*xa*xc)-xa*(3*ac*(xa+xb)+xa*cb))/4;
        xy = xa*ya*(6*ab+6*ac-3*bc+cb); ac = ya*ya; cb = xa*xa;
        xy = 3*(xy+9*f*(cb*yb*yc-xb*xc*ac)-18*xb*yb*ab)/8;

        if (ex < 0)
        {
            dx = -dx; dy = -dy; xx = -xx; yy = -yy; xy = -xy; ac = -ac; cb = -cb;
        }

        ab = 6*ya*ac; ac = -6*xa*ac; bc = 6*ya*cb; cb = -6*xa*cb;
        dx += xy; ex = dx+dy; dy += xy;

        for (pxy = &xy, fx = fy = f; x0 != x3 && y0 != y3; )
        {
            gc_draw_px(x0,y0);
            do {
                if (dx > *pxy || dy < *pxy) goto exit;
                y1 = 2*ex-dy;
                if (2*ex >= dx) {
                    fx--; ex += dx += xx; dy += xy += ac; yy += bc; xx += ab;
                }
                if (y1 <= 0) {
                    fy--; ex += dy += yy; dx += xy += bc; xx += ac; yy += cb;
                }
            } while (fx > 0 && fy > 0);
            if (2*fx <= f) { x0 += sx; fx += f; }
            if (2*fy <= f) { y0 += sy; fy += f; }
            if (pxy == &xy && dx < 0 && dy > 0) pxy = &EP;
        }
        exit: xx = x0; x0 = x3; x3 = xx; sx = -sx; xb = -xb;
        yy = y0; y0 = y3; y3 = yy; sy = -sy; yb = -yb; x1 = x2;
    } while (leg--);
    gc_draw_line(x0,y0, x3,y3);
}

void gc_draw_cubic_bezier (int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3)
{
    int n = 0, i = 0;
    long xc = x0+x1-x2-x3, xa = xc-4*(x1-x2);
    long xb = x0-x1-x2+x3, xd = xb+4*(x1+x2);
    long yc = y0+y1-y2-y3, ya = yc-4*(y1-y2);
    long yb = y0-y1-y2+y3, yd = yb+4*(y1+y2);
    float fx0 = x0, fx1, fx2, fx3, fy0 = y0, fy1, fy2, fy3;
    double t1 = xb*xb-xa*xc, t2, t[5];

    if (xa == 0) {
        if (abs(xc) < 2*abs(xb)) t[n++] = xc/(2.0*xb);
    } else if (t1 > 0.0) {
        t2 = sqrt(t1);
        t1 = (xb-t2)/xa; if (fabs(t1) < 1.0) t[n++] = t1;
        t1 = (xb+t2)/xa; if (fabs(t1) < 1.0) t[n++] = t1;
    }

    t1 = yb*yb-ya*yc;

    if (ya == 0) {
        if (abs(yc) < 2*abs(yb)) t[n++] = yc/(2.0*yb);
    } else if (t1 > 0.0) {
        t2 = sqrt(t1);
        t1 = (yb-t2)/ya; if (fabs(t1) < 1.0) t[n++] = t1;
        t1 = (yb+t2)/ya; if (fabs(t1) < 1.0) t[n++] = t1;
    }

    for (i = 1; i < n; i++) if ((t1 = t[i-1]) > t[i]) { t[i-1] = t[i]; t[i] = t1; i = 0; }
    t1 = -1.0; t[n] = 1.0;

    for (i = 0; i <= n; i++)
    {
        t2 = t[i];
        fx1 = (t1*(t1*xb-2*xc)-t2*(t1*(t1*xa-2*xb)+xc)+xd)/8-fx0;
        fy1 = (t1*(t1*yb-2*yc)-t2*(t1*(t1*ya-2*yb)+yc)+yd)/8-fy0;
        fx2 = (t2*(t2*xb-2*xc)-t1*(t2*(t2*xa-2*xb)+xc)+xd)/8-fx0;
        fy2 = (t2*(t2*yb-2*yc)-t1*(t2*(t2*ya-2*yb)+yc)+yd)/8-fy0;
        fx0 -= fx3 = (t2*(t2*(3*xb-t2*xa)-3*xc)+xd)/8;
        fy0 -= fy3 = (t2*(t2*(3*yb-t2*ya)-3*yc)+yd)/8;
        x3 = floor(fx3+0.5); y3 = floor(fy3+0.5);
        if (fx0 != 0.0) { fx1 *= fx0 = (x0-x3)/fx0; fx2 *= fx0; }
        if (fy0 != 0.0) { fy1 *= fy0 = (y0-y3)/fy0; fy2 *= fy0; }
        if (x0 != x3 || y0 != y3) cubicBezierSeg(x0,y0, x0+fx1,y0+fy1, x0+fx2,y0+fy2, x3,y3);
        x0 = x3; y0 = y3; fx0 = fx3; fy0 = fy3; t1 = t2;
    }
}

#define MAX_STACK_SIZE 5000
typedef struct {
    int x;
    int y;
} Point;
Point* colorStack;
int top = -1;

void pushToColorStack (int x, int y)
{
    if (top < MAX_STACK_SIZE - 1)
    {
        Point p;
        p.x = x;
        p.y = y;
        colorStack[++top] = p;
    }
}

Point popFromColorStack ()
{
    Point p;
    p.x = -1;
    p.y = -1;

    if (top >= 0) p = colorStack[top--];
    return p;
}

int isEmpty ()
{
    return top == -1;
}

static bool colorsMatch (gc_color_t color1, gc_color_t color2)
{
    if (color1.r == color2.r && color1.g == color2.g && color1.b == color2.b) return true;
    return false;
}

static bool isValidPixelPath (int x, int y)
{
    if (x >= 0 && x <= GC_DISPLAY_WIDTH && y >= 0 && y <= GC_DISPLAY_HEIGHT)
    {
        if (!colorsMatch(currentColor, buffer_get_px(x, y))) return true;
    }

    return false;
}

void gc_draw_fill (int startX, int startY)
{
    colorStack = heap_caps_malloc(sizeof(Point) * MAX_STACK_SIZE, MALLOC_CAP_SPIRAM);
    pushToColorStack(startX, startY);
    while (!isEmpty())
    {
        Point current = popFromColorStack();
        int x = current.x;
        int y = current.y;

        if (colorsMatch(currentColor, buffer_get_px(x, y))) continue;
        gc_draw_px(x, y);

        if (isValidPixelPath(x - 1, y)) pushToColorStack(x - 1, y);
        if (isValidPixelPath(x + 1, y)) pushToColorStack(x + 1, y);
        if (isValidPixelPath(x, y - 1)) pushToColorStack(x, y - 1);
        if (isValidPixelPath(x, y + 1)) pushToColorStack(x, y + 1);
    }
    free(colorStack);
}

void gc_display_clear ()
{
    for (int x = 0; x < GC_DISPLAY_WIDTH; x++)
    {
        for (int y = 0; y < GC_DISPLAY_HEIGHT; y++)
        {
            gc_draw_px(x, y);
        }
    }
}