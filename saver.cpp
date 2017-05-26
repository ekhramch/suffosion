#include "saver.hpp"

//---------------------------------------------------------------------------
template <typename T>
std::string to_string(T v)
{
    std::ostringstream ss;//create a stringstream
    ss << v;//add number to the stream
    return ss.str();//return a string with the contents of the stream    
}

//---------------------------------------------------------------------------
saver::saver(const std::string &fname,
        hsize_t nx, hsize_t ny, hsize_t nz,
        double h
        )
    : hdf(fname + ".h5", H5F_ACC_TRUNC), dtype(H5::PredType::NATIVE_DOUBLE),
      xmf_name(fname + ".xmf"), step(0), nx(nx), ny(ny), nz(nz), hdf_name(fname + ".h5")
{
    std::vector<double> x(nx), y(ny), z(nz);
    for(int i = 0; i < nx; ++i) x[i] = i * h;
    for(int i = 0; i < ny; ++i) y[i] = i * h;
    for(int i = 0; i < nz; ++i) z[i] = i * h;

    {
        H5::DataSpace space(1, &nx);
        H5::DataSet ds = hdf.createDataSet("x", dtype, space);
        ds.write(x.data(), dtype);
    }
    {
        H5::DataSpace space(1, &ny);
        H5::DataSet ds = hdf.createDataSet("y", dtype, space);
        ds.write(y.data(), dtype);
    }
    {
        H5::DataSpace space(1, &nz);
        H5::DataSet ds = hdf.createDataSet("z", dtype, space);
        ds.write(z.data(), dtype);
    }

    {
        hsize_t zero  = 0;
        hsize_t unlim = H5S_UNLIMITED;
        H5::DataSpace space(1, &zero, &unlim);

        hsize_t chunk_size = 1024;
        H5::DSetCreatPropList p;
        p.setChunk(1, &chunk_size);

        time_ds = hdf.createDataSet("time", dtype, space, p);
    }

    {
        hsize_t dim[] = {nz, ny, nx};
        fspace = H5::DataSpace(3, dim);
    }

    hsize_t chunk_size[] = {
        std::min<uint>(32, nz),
        std::min<uint>(32, ny),
        std::min<uint>(32, nx)
    };
    prop.setChunk(3, chunk_size);
    prop.setDeflate(9);

    xmf.append_child(pugi::node_doctype).set_value("Xdmf SYSTEM \"Sdmf.dtd\" []");

    pugi::xml_node xdmf = xmf.append_child("Xdmf");
    xdmf.append_attribute("xmlns:xi").set_value("http://www.w3.org/2003/XInclude");
    xdmf.append_attribute("Version" ).set_value("2.2");
    pugi::xml_node domain = xdmf.append_child("Domain");

    time_collection = domain.append_child("Grid");
    time_collection.append_attribute("Name"          ).set_value("Time");
    time_collection.append_attribute("GridType"      ).set_value("Collection");
    time_collection.append_attribute("CollectionType").set_value("Temporal");
}

//---------------------------------------------------------------------------
void saver::add_step(double time, const std::map<std::string, double*> &v) {
    ++step;

    {
        hsize_t cur_size;
        time_ds.getSpace().getSimpleExtentDims(&cur_size);

        if (step >= cur_size)
            time_ds.extend(&step);

        hsize_t start = step - 1;
        hsize_t count = 1;

        H5::DataSpace mspace = H5::DataSpace(1, &count);
        H5::DataSpace fspace = time_ds.getSpace();

        fspace.selectHyperslab(H5S_SELECT_SET, &count, &start);

        time_ds.write(&time, dtype, mspace, fspace);
    } 

    pugi::xml_node grid = time_collection.append_child("Grid");
    grid.append_child("Time").append_attribute("Value").set_value(to_string(time).c_str());
    grid.append_attribute("Name"    ).set_value("data");
    grid.append_attribute("GridType").set_value("Uniform");

    pugi::xml_node topo = grid.append_child("Topology");
    topo.append_attribute("TopologyType"    ).set_value("3DRectMesh");
    topo.append_attribute("NumberOfElements").set_value( (
                to_string(nz) + " " +
                to_string(ny) + " " +
                to_string(nx)
                ).c_str());

    pugi::xml_node geom = grid.append_child("Geometry");
    geom.append_attribute("GeometryType").set_value("VXVYVZ");

    {
        pugi::xml_node data = geom.append_child("DataItem");
        data.append_attribute("Dimensions").set_value(to_string(nx).c_str());
        data.append_attribute("NumberType").set_value("Float");
        data.append_attribute("Precision" ).set_value("8");
        data.append_attribute("Format"    ).set_value("HDF");

        data.append_child(pugi::node_pcdata).set_value((hdf_name + ":/x").c_str());
    }

    {
        pugi::xml_node data = geom.append_child("DataItem");
        data.append_attribute("Dimensions").set_value(to_string(ny).c_str());
        data.append_attribute("NumberType").set_value("Float");
        data.append_attribute("Precision" ).set_value("8");
        data.append_attribute("Format"    ).set_value("HDF");

        data.append_child(pugi::node_pcdata).set_value((hdf_name + ":/y").c_str());
    }

    {
        pugi::xml_node data = geom.append_child("DataItem");
        data.append_attribute("Dimensions").set_value(to_string(nz).c_str());
        data.append_attribute("NumberType").set_value("Float");
        data.append_attribute("Precision" ).set_value("8");
        data.append_attribute("Format"    ).set_value("HDF");

        data.append_child(pugi::node_pcdata).set_value((hdf_name + ":/z").c_str());
    }

    std::ostringstream group;
    group << "step-" << step;
    hdf.createGroup(group.str());

    for(auto i : v) {
        std::ostringstream dsname;
        dsname << "step-" << step << "/" << i.first;
        hdf.createDataSet(dsname.str(), dtype, fspace, prop).write(i.second, dtype);


        pugi::xml_node attr = grid.append_child("Attribute");

        attr.append_attribute("Name"         ).set_value(i.first.c_str());
        attr.append_attribute("AttributeType").set_value("Scalar");
        attr.append_attribute("Center"       ).set_value("Node");

        pugi::xml_node data = attr.append_child("DataItem");
        data.append_attribute("Dimensions").set_value( (
                    to_string(nz) + " " +
                    to_string(ny) + " " +
                    to_string(nx)
                    ).c_str());
        data.append_attribute("NumberType").set_value("Float");
        data.append_attribute("Precision" ).set_value("8");
        data.append_attribute("Format"    ).set_value("HDF");

        data.append_child(pugi::node_pcdata).set_value( (
                    hdf_name + ":/step-" + to_string(step) + "/" + i.first
                    ).c_str());
    }

    xmf.save_file(xmf_name.c_str(), "  ");
}
