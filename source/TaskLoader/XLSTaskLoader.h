#ifndef XLSTASKLOADER_H
#define XLSTASKLOADER_H

#include "../../include/TaskLoader.h"

class XLSTaskLoader : public ITaskLoader {
public:
    XLSTaskLoader(const char* _pathToTaskFile) : ITaskLoader(_pathToTaskFile) { }

    int load(std::size_t numberOfPoints);
};

#endif // XLSTASKLOADER_H
