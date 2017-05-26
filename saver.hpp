#ifndef SAVER_HPP
#define SAVER_HPP

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <H5Cpp.h>
#include <pugixml.hpp>

class saver {
    public:
        saver(const std::string &fname, hsize_t nx, hsize_t ny, hsize_t nz, double h);
        void add_step(double time, const std::map<std::string, double*> &v);

    private:
        H5::H5File hdf;
        H5::DataSet time_ds;
        H5::DataSpace fspace;
        H5::DSetCreatPropList prop;


        std::string xmf_name;
        pugi::xml_document xmf;
        pugi::xml_node time_collection;
        H5::FloatType dtype;
        hsize_t step;
        hsize_t nx;
        hsize_t ny;
        hsize_t nz;
        std::string hdf_name;
};

#endif
