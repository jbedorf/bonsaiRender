#ifndef _RENDERLOOP_H_
#define _RENDERLOOP_H_

void initGL(int argc, char** argv, const char *fullScreenMode, bool &stereo);
#if 0
void initAppRenderer(int argc, char** argv, void *tree, 
                     octree::IterationData &idata,
                     bool showFPS, bool stereo);
#endif
#endif // _RENDERLOOP_H_
