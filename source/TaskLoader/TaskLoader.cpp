#include "../../include/TaskLoader.h"
#include "XLSTaskLoader.h"
#include "TXTTaskLoader.h"

#include "../../lib/libxl/libxl.h"

#include <fstream>

using libxl::Book;
using libxl::Sheet;

ITaskLoader::ITaskLoader(const char* _pathToTaskFile) : pathToTaskFile(_pathToTaskFile) { }
std::vector<double>& ITaskLoader::getYears()   { return years;   }
std::vector<double>& ITaskLoader::getHeights() { return heights; }

int TXTTaskLoader::load(size_t numberOfPoints) {
    years.resize(numberOfPoints);
    heights.resize(numberOfPoints);

    std::ifstream taskFile;
    taskFile.open(pathToTaskFile, std::ifstream::in);
    if(!taskFile.is_open()) {
        return TXT_NOT_LOADED;
    }

    for(int i = 0; i < numberOfPoints; i++) {
        taskFile >> years[i] >> heights[i];
    }

    taskFile.close();

    return LOAD_OK;
}

int XLSTaskLoader::load(size_t numberOfPoints) {
    years.reserve(numberOfPoints);
    heights.reserve(numberOfPoints);

    Book* book = xlCreateBook();
    book->load(pathToTaskFile);
    if(!book) {
        return XLS_BOOK_NOT_LOADED;
    }

    Sheet* sheet = book->getSheet(0);
    if(!sheet) {
        return XLS_SHEET_NOT_LOADED;
    }

    for(int i = 0; i < numberOfPoints; i++) {
        years.push_back(sheet->readNum(i, 0));
        heights.push_back(sheet->readNum(i, 1));
    }

    return LOAD_OK;
}

ITaskLoader* TaskLoaderFactory::makeLoader(TaskFileFormat format, const char* pathToTaskFile) {
    switch (format) {
        case XLSLib: return new XLSTaskLoader(pathToTaskFile); break;
        case TXT: return new TXTTaskLoader(pathToTaskFile); break;
        default: break;
    }
}
