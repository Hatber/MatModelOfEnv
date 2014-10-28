#ifndef TXTFILELOADER_H
#define TXTFILELOADER_H

#include "../../include/TaskLoader.h"

class TXTTaskLoader : public ITaskLoader {
public:
    TXTTaskLoader(const char* _pathToTaskFile) : ITaskLoader(_pathToTaskFile) { }

    int load(std::size_t numberOfPoints);
};

#endif // TXTFILELOADER_H
