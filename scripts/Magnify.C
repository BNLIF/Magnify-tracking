void Magnify(const char* filename="../data/track_com_mc.root", int sign=0);
void Magnify(const char* filename, int sign)
{
    // Data *data = new Data(filename);

    GuiController *gc = new GuiController(
        gClient->GetRoot(),
        1600,
        900,
        filename,
        sign
    );


}
