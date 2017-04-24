/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataMan.cpp
 *
 *  Created on: Apr 12, 2017
 *      Author: Jason Wang
 */

#include "DataMan.h"

int DataMan::init(json a_jmsg) { return 0; }

int DataMan::put(const void *a_data, std::string p_doid, std::string p_var,
                 std::string p_dtype, std::vector<size_t> p_putshape,
                 std::vector<size_t> p_varshape, std::vector<size_t> p_offset,
                 size_t p_timestep, int p_tolerance, int p_priority)
{
    return DataMan::put(a_data, p_doid, p_var, p_dtype, p_putshape, p_varshape,
                        p_offset, p_timestep, p_tolerance, p_priority);
}

int DataMan::put(const void *a_data, json a_jmsg)
{
    if (m_cache_size > 0)
    {
        // to be implemented
    }
    else
    {
        put_begin(a_data, a_jmsg);
        put_end(a_data, a_jmsg);
    }

    return 0;
}

void DataMan::add_file(std::string p_method) {}

void DataMan::add_stream(json a_jmsg)
{

    std::string method;

    if (a_jmsg["method"].is_string())
    {
        method = a_jmsg["method"];
    }

    logging("Streaming method " + method + " added");

    if (a_jmsg["cachesize"].is_number())
    {
        m_cache_size = a_jmsg["cachesize"].get<size_t>();
    }

    if (m_tolerance.size() < m_num_channels)
    {
        for (int i = 0; i < m_num_channels; i++)
        {
            m_tolerance.push_back(0);
        }
    }
    if (m_priority.size() < m_num_channels)
    {
        for (int i = 0; i < m_num_channels; i++)
        {
            m_priority.push_back(100 / (i + 1));
        }
    }

    auto man = get_man(method);
    if (man)
    {
        man->init(a_jmsg);
        this->add_next(method, man);
    }
    if (a_jmsg["compression_method"].is_string())
    {
        if (a_jmsg["compression_method"] != "null")
        {
            add_man_to_path(a_jmsg["compression_method"], method, a_jmsg);
        }
    }
}

void DataMan::flush() { flush_next(); }

int DataMan::get(void *a_data, json &a_jmsg) { return 0; }