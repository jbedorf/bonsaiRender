#ifndef _RENDERLOOP_H_
#define _RENDERLOOP_H_

#include "RendererData.h"

void initGL(int argc, char** argv, const char *fullScreenMode, bool &stereo);
void initAppRenderer(int argc, char** argv, 
                     RendererData &data,
                     bool showFPS, bool stereo);
#endif // _RENDERLOOP_H_
