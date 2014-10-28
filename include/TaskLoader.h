#ifndef ITASKLOADER_H
#define ITASKLOADER_H

#include <vector>

const int LOAD_OK = 0;
const int TXT_NOT_LOADED = -1;
const int XLS_BOOK_NOT_LOADED = -2;
const int XLS_SHEET_NOT_LOADED = -3;

enum TaskFileFormat{XLSLib, TXT};

class ITaskLoader {
public:
    ITaskLoader(const char* _pathToTaskFile);
    virtual ~ITaskLoader() { }

    virtual int load(std::size_t numberOfPoints) = 0;

    std::vector<double>& getYears();
    std::vector<double>& getHeights();
protected:
    const char* pathToTaskFile;

    std::vector<double> years;
    std::vector<double> heights;
};

class TaskLoaderFactory {
    public : static ITaskLoader* makeLoader(TaskFileFormat format, const char* pathToTaskFile);
};

#endif // ITASKLOADER_H
