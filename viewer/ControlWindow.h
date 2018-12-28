#ifndef CONTROL_WINDOW_H
#define CONTROL_WINDOW_H

#include "TGFrame.h"

class TGNumberEntry;
class TGCheckButton;
class TGTextButton;

class ControlWindow: public TGHorizontalFrame
{
public:
    ControlWindow(const TGWindow *p, int w, int h);
    virtual ~ControlWindow();

    TGNumberEntry *clusterEntry;
    TGNumberEntry *zoomEntry;
    TGNumberEntry *clusterIdEntry;


    // TGNumberEntry *threshEntry[3];
    // TGNumberEntry *zAxisRangeEntry[2];
    // TGNumberEntry *timeRangeEntry[2];
    // TGNumberEntry *adcRangeEntry[2];
    TGCheckButton *badChanelButton;
    TGTextButton  *unZoomButton;

    // TGTextButton  *setThreshButton;
    // TGNumberEntry *threshScaleEntry;

    ClassDef(ControlWindow, 0)
};

#endif
