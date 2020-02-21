#include "ControlWindow.h"

#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGButtonGroup.h"
#include "TGListBox.h"
#include "TGLabel.h"
// #include "TGDoubleSlider.h"
#include "TString.h"

#include <iostream>
using namespace std;

ControlWindow::ControlWindow(const TGWindow *p, int w, int h)
    :TGHorizontalFrame(p, w, h)
{
    TGGroupFrame *group_general = new TGGroupFrame(this, "General", kHorizontalFrame);
    group_general->SetTitlePos(TGGroupFrame::kLeft);
    AddFrame(group_general, new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 1, 1));

    group_general->AddFrame(new TGLabel(group_general, "cluster index: "), new TGLayoutHints(kLHintsTop | kLHintsLeft,  2, 2, 1, 1));
    clusterEntry = new TGNumberEntry(group_general, 0, 5, -1,
        TGNumberFormat::kNESInteger,
        TGNumberFormat::kNEANonNegative,
        TGNumberFormat::kNELLimitMinMax,
        0, 20);
    group_general->AddFrame(clusterEntry, new TGLayoutHints(kLHintsTop | kLHintsLeft,  1, 1, 1, 1));

    group_general->AddFrame(new TGLabel(group_general, "cluster id: "), new TGLayoutHints(kLHintsTop | kLHintsLeft,  2, 2, 1, 1));
    clusterIdEntry = new TGNumberEntry(group_general, 0, 5, -1,
        TGNumberFormat::kNESInteger,
        TGNumberFormat::kNEANonNegative,
        TGNumberFormat::kNELLimitMinMax,
        0, 1000);
    group_general->AddFrame(clusterIdEntry, new TGLayoutHints(kLHintsTop | kLHintsLeft,  1, 1, 1, 1));

    badChanelButton = new TGCheckButton(group_general, "bad channel  ");
    badChanelButton->SetState(kButtonUp);
    group_general->AddFrame(badChanelButton, new TGLayoutHints(kLHintsTop | kLHintsLeft,  1, 1, 1, 1));

    drawTrackButton = new TGCheckButton(group_general, "draw track  ");
    drawTrackButton->SetState(kButtonDown);
    group_general->AddFrame(drawTrackButton, new TGLayoutHints(kLHintsTop | kLHintsLeft,  1, 1, 1, 1));


    TGGroupFrame *group_misc = new TGGroupFrame(this, "Zoom", kHorizontalFrame);
    group_misc->SetTitlePos(TGGroupFrame::kLeft);
    AddFrame(group_misc, new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 1, 1));

    group_misc->AddFrame(new TGLabel(group_misc, "point index: "), new TGLayoutHints(kLHintsTop | kLHintsLeft,  2, 2, 1, 1));
    pointIndexEntry = new TGNumberEntry(group_misc, 2, 5, -1,
        TGNumberFormat::kNESInteger,
        TGNumberFormat::kNEANonNegative,
        TGNumberFormat::kNELLimitMinMax,
        1, 1000);
    group_misc->AddFrame(pointIndexEntry, new TGLayoutHints(kLHintsTop | kLHintsLeft,  1, 1, 1, 1));

    group_misc->AddFrame(new TGLabel(group_misc, "zoom +- 10x: "), new TGLayoutHints(kLHintsTop | kLHintsLeft,  2, 2, 1, 1));
    zoomEntry = new TGNumberEntry(group_misc, 2, 5, -1,
        TGNumberFormat::kNESInteger,
        TGNumberFormat::kNEANonNegative,
        TGNumberFormat::kNELLimitMinMax,
        1, 100);
    group_misc->AddFrame(zoomEntry, new TGLayoutHints(kLHintsTop | kLHintsLeft,  1, 1, 1, 1));


    unZoomButton = new TGTextButton(group_misc, "UnZoom");
    group_misc->AddFrame(unZoomButton, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 3, 2, 2, 2));


    TGGroupFrame *group_range = new TGGroupFrame(this, "Range", kHorizontalFrame);
    group_range->SetTitlePos(TGGroupFrame::kLeft);
    AddFrame(group_range, new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 1, 1));

    TString str[4] = {"time: ", "u: ", "v: ", "w: "};
    int min[4] = {0, 0, 2400, 4800};
    int max[4] = {2400, 2400, 4800, 8256};

    for (int i=0; i<4; i++) {
        group_range->AddFrame(new TGLabel(group_range, str[i].Data()), new TGLayoutHints(kLHintsTop | kLHintsLeft,  2, 2, 1, 1));

        minEntry[i] = new TGNumberEntry(group_range, min[i], 5, -1,
            TGNumberFormat::kNESInteger,
            TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax,
            min[i], max[i]);
        group_range->AddFrame(minEntry[i], new TGLayoutHints(kLHintsTop | kLHintsLeft,  1, 1, 1, 1));

        maxEntry[i] = new TGNumberEntry(group_range, max[i], 5, -1,
            TGNumberFormat::kNESInteger,
            TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax,
            min[i], max[i]);
        group_range->AddFrame(maxEntry[i], new TGLayoutHints(kLHintsTop | kLHintsLeft,  1, 1, 1, 1));
    }

    rangeZoomButton = new TGTextButton(group_range, "Zoom");
    group_range->AddFrame(rangeZoomButton, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 3, 2, 2, 2));

    keepRangeButton = new TGCheckButton(group_range, "keep");
    keepRangeButton->SetState(kButtonUp);
    group_range->AddFrame(keepRangeButton, new TGLayoutHints(kLHintsTop | kLHintsLeft,  1, 1, 1, 1));

}

ControlWindow::~ControlWindow()
{
}
