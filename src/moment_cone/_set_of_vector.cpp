#include <vector>
#include <set>
#include <cstdint>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

template <typename T>
class SetOfVector
{
private:
    class comparison;

public:
    class const_iterator;

public:
    using value_type = const typename std::vector<T>;

    SetOfVector(std::size_t width)
        : set(comparison(this))
        , width(width)
    {}

    bool add(value_type const& element)
    {
        for (T const& value : element) // size should be checked before in Python
            data.push_back(value);

        auto [it, added] = set.insert(set.size());
        if (not added)
            data.resize(data.size() - width);
        return added;
    }

    std::size_t size() const { return set.size(); }
    const_iterator begin() const { return const_iterator(this); }
    const_iterator end() const { return const_iterator(this, size()); }

    void clear()
    {
        set.clear();
        data.clear();
    }

public:
    std::set<std::size_t, comparison> set;
    std::vector<T> data;
    std::size_t width;
};

template <typename T>
class SetOfVector<T>::comparison
{
private:
    SetOfVector<T> const* sov;

public:
    explicit comparison(SetOfVector<T> const* sov) : sov(sov) {}
    bool operator() (std::size_t i, std::size_t j) const
    {
        T const* lhs = sov->data.data() + i * sov->width;
        T const* rhs = sov->data.data() + j * sov->width;
        for (std::size_t pos = 0; pos < sov->width; ++pos)
        {
            if (lhs[pos] != rhs[pos])
                return lhs[pos] < rhs[pos];
        }
        return false;
    }
};

template <typename T>
class SetOfVector<T>::const_iterator
{
private:
    SetOfVector<T> const* sov;
    std::size_t idx;

public:
    using value_type = typename SetOfVector<T>::value_type;
    using distance_type = std::ptrdiff_t;
    using pointer = const value_type*;
    using reference = value_type;

    explicit const_iterator(SetOfVector<T> const* sov, std::size_t idx = 0)
        : sov(sov), idx(idx) {}
    bool operator== (const_iterator const& other) const { return idx == other.idx; }
    bool operator!= (const_iterator const& other) const { return idx != other.idx; }
    const_iterator& operator++ () { ++idx; return *this; }
    const_iterator operator++ (int) { const_iterator copy(*this); ++(*this); return copy; }
    reference operator* () const
    {
        return reference(
            sov->data.begin() + idx * sov->width,
            sov->data.begin() + (idx + 1) * sov->width
        );
    }
};


#define AddSetOfVector(name, type) \
    py::class_<SetOfVector<type>>(m, name) \
        .def(py::init<std::size_t>()) \
        .def("add", &SetOfVector<type>::add) \
        .def("clear", &SetOfVector<type>::clear) \
        .def("__len__", &SetOfVector<type>::size) \
        .def("__iter__", [](SetOfVector<type> const& sov) { return py::make_iterator(sov.begin(), sov.end()); }, py::keep_alive<0, 1>());

PYBIND11_MODULE(_set_of_vector, m) {
    AddSetOfVector("SetOfVector64", std::int_least64_t)
    AddSetOfVector("SetOfVector32", std::int_least32_t)
    AddSetOfVector("SetOfVector16", std::int_least16_t)
    AddSetOfVector("SetOfVector8", std::int_least8_t)
}
