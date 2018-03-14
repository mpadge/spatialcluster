#pragma once

struct edge_component
{
    // used only for edge sampling on graphs without component numbers
    edge_id_t edge;
    unsigned int component;
};

struct vertex_t
{
    private:
        std::unordered_set <vertex_id_t> in, out;

    public:
        void add_neighbour_in (vertex_id_t vert_id) { in.insert (vert_id); }
        void add_neighbour_out (vertex_id_t vert_id) { out.insert (vert_id); }
        unsigned long int get_degree_in () { return in.size (); }
        unsigned long int get_degree_out () { return out.size (); }

        std::unordered_set <vertex_id_t> get_all_neighbours ()
        {
            std::unordered_set <vertex_id_t> all_neighbours = in;
            all_neighbours.insert (out.begin (), out.end ());
            return all_neighbours;
        }

        std::unordered_set <vertex_id_t> get_in_neighbours ()
        {
            return in;
        }

        std::unordered_set <vertex_id_t> get_out_neighbours ()
        {
            return out;
        }
};

struct edge_t
{
    private:
        vertex_id_t from, to;
        edge_id_t id;

    public:
        vertex_id_t get_from_vertex () { return from; }
        vertex_id_t get_to_vertex () { return to; }
        edge_id_t getID () { return id; }

        edge_t (vertex_id_t from_id, vertex_id_t to_id, edge_id_t id_in)
        {
            this -> to = to_id;
            this -> from = from_id;
            this -> id = id_in;
        }
};


typedef std::unordered_map <vertex_id_t, vertex_t> vertex_map_t;
typedef std::unordered_map <edge_id_t, edge_t> edge_map_t;
typedef std::unordered_map <vertex_id_t,
        std::unordered_set <edge_id_t> > vert2edge_map_t;

void add_to_v2e_map (vert2edge_map_t &vert2edge_map, const vertex_id_t vid,
        const edge_id_t eid);

void erase_from_v2e_map (vert2edge_map_t &vert2edge_map, const vertex_id_t vid,
        const edge_id_t eid);

bool graph_has_components (const Rcpp::DataFrame &graph);

void graph_from_df (const Rcpp::DataFrame &gr, vertex_map_t &vm,
        edge_map_t &edge_map, vert2edge_map_t &vert2edge_map);

unsigned int identify_graph_components (vertex_map_t &v,
        std::unordered_map <vertex_id_t, unsigned int> &com);

Rcpp::List rcpp_get_component_vector (const Rcpp::DataFrame &graph);
