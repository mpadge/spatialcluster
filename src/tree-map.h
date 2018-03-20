#pragma once

#include <map>
#include <algorithm> // min_element
#include <iostream>

template <typename T1, typename T2>
class TreeMap {
    typedef std::map <T1, T2> TreeMap_t;
    typedef std::pair <T1, T2> TreePair_t;
    private:
        TreeMap_t tree_map;

    public:
        void tree_insert (T1 key, T2 value)
        {
            tree_map.emplace (key, value);
        }

        void tree_delete (T1 key)
        {
            tree_map.erase (key);
        }

        static bool compare (TreePair_t i, TreePair_t j) 
        { 
            return i.second < j.second; 
        }

        double getMin ()
        {
            TreePair_t min = *min_element (tree_map.begin(), tree_map.end(),
                    compare);
            return min.second; 
        }
};
