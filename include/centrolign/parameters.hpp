#ifndef centrolign_parameters_hpp
#define centrolign_parameters_hpp

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <stdexcept>
#include <type_traits>

#include "centrolign/logging.hpp"
#include "centrolign/core.hpp"

namespace centrolign {

/*
 * The command line parameters and their defaults
 */
struct Parameters {
public:
    Parameters();
    Parameters(std::istream& in);
    ~Parameters() = default;
    
    // generate a config YAML that would recreate this set of parameters
    std::string generate_config() const;
    
    // check logical validity of parameters, throw an error if anything is off
    void validate() const;
    
    // pass all of the parameters to configurable modules
    void apply(Core& core) const;
    
    // set the parameter of this name to the given value
    template<typename T>
    void set(const std::string& param_name, T value);
    
    template<typename T>
    T get(const std::string& param_name) const;
    
    bool operator==(const Parameters& other) const;
    bool operator!=(const Parameters& other) const;
    
private:
    
    // functional submodules of the algorithm
    enum submodule_t {
        IO,
        MatchFinding,
        Anchoring,
        IdentifyingAlignability,
        Aligning,
        InducingCycles,
        DeveloperTools
    };
    
    // data types
    enum type_t {
        Integer,
        Bool,
        Double,
        Enum,
        String,
        DoubleArray3,
        IntegerArray3
    };
    
    // a union that could represent any of the data types
    union value_t {
        int64_t i;
        bool b;
        double d;
        std::string* s;
        std::array<double, 3> da;
        std::array<int64_t, 3> ia;
    };
    
    /*
     * A single parameter
     */
    struct Parameter {
    public:
        
        template<typename T>
        Parameter(type_t type, const std::string& name, const std::string& help, T value);
        Parameter() = default;
        Parameter(Parameter&& other);
        Parameter(const Parameter& other);
        ~Parameter();
        
        Parameter& operator=(Parameter&& other);
        Parameter& operator=(const Parameter& other);
        
        template<typename T>
        void set(T value);
        
        template<typename T>
        T get() const;
        
        type_t get_type() const;
        const std::string& get_name() const;
        const std::string& get_help() const;
        std::string value_str() const;
        
    private:
        
        template<typename T>
        void check_type() const;
        
        template<typename T>
        void set_internal(T val);
        
        template<typename T>
        T get_internal() const;
        
        const type_t type = (type_t) -1;
        const std::string name;
        const std::string help;
        value_t value = { 0 };
        
    };
    
    // the parameters, organized by their submodule
    std::map<submodule_t, std::pair<std::string, std::vector<Parameter>>> params;
    // the location of the parameter with this name within the main parameter container
    std::unordered_map<std::string, std::pair<submodule_t, size_t>> param_position;
    
    // set up methods
    void initialize_submodule(submodule_t module, const std::string& description);
    template<typename T>
    void add_parameter(submodule_t module, const std::string& name, type_t type, T default_value, const std::string& help);
    
    // access a parameter by name
    const Parameter& parameter(const std::string& name) const;
    Parameter& parameter(const std::string& name);
    
    template<typename T>
    void enforce_lim(const std::string& name, bool lim_is_min, bool equal_ok, T lim) const;
    template<typename T>
    void enforce_lt(const std::string& name, T upper_bound) const;
    template<typename T>
    void enforce_leq(const std::string& name, T upper_bound) const;
    template<typename T>
    void enforce_gt(const std::string& name, T lower_bound) const;
    template<typename T>
    void enforce_geq(const std::string& name, T lower_bound) const;
    template<typename T>
    void enforce_range(const std::string& name, T lower_bound, T upper_bound) const;
};






/*
 * Template implementations
 */

template<typename T>
void Parameters::add_parameter(submodule_t module, const std::string& name, type_t type, T default_value, const std::string& help) {
    
    if (param_position.count(name)) {
        throw std::runtime_error("Added parameter with duplicate name " + name);
    }
    if (!params.count(module)) {
        throw std::runtime_error("Added parameter to nonexistent submodule " + std::to_string((int) module));
    }
    param_position[name] = std::make_pair(module, params.at(module).second.size());
    params[module].second.emplace_back(type, name, help, default_value);
}

template<typename T>
void Parameters::set(const std::string& param_name, T value) {
    parameter(param_name).set<T>(value);
}

template<typename T>
T Parameters::get(const std::string& param_name) const {
    return parameter(param_name).get<T>();
}

template<typename T>
Parameters::Parameter::Parameter(type_t type, const std::string& name, const std::string& help, T value) : type(type), name(name), help(help) {
    set<T>(value);
}

template<typename T>
void Parameters::Parameter::check_type() const {
    static const std::map<type_t, std::string> type_names = {
        {Integer, "integer"},
        {Double, "double"},
        {Bool, "bool"},
        {Enum, "enum"},
        {String, "string"},
        {DoubleArray3, "double[3]"},
        {IntegerArray3, "integer[3]"}
    };
    if ((type == Integer && !std::is_integral<T>::value) ||
        (type == Double && !std::is_floating_point<T>::value) ||
        (type == Bool && !std::is_integral<T>::value) ||
        (type == Enum && !std::is_enum<T>::value && !std::is_integral<T>::value) || // hard not to slip integers on the setting side
        (type == String && !std::is_same<T, const char*>::value  && !std::is_same<T, char*>::value && !std::is_same<T, std::string>::value) ||
        (type == DoubleArray3 && std::is_array<T>::value) ||
        (type == IntegerArray3 && std::is_array<T>::value)) {
        throw std::runtime_error(std::string("Value of type ") + typeid(T).name() + " for parameter " + name + " does not match expected type " + type_names.at(type));
    }
}

template<typename T>
void Parameters::Parameter::set(T value) {
    check_type<T>();
    set_internal<T>(value);
}

// default
template<typename T>
inline void Parameters::Parameter::set_internal(T value) {
    throw std::runtime_error(std::string("Parameter setting not implemented for type ") + typeid(value).name());
}
// specializations
template<>
inline void Parameters::Parameter::set_internal<int64_t>(int64_t value) {
    this->value.i = value;
}
template<>
inline void Parameters::Parameter::set_internal<int32_t>(int32_t value) {
    this->value.i = value;
}
template<>
inline void Parameters::Parameter::set_internal<uint64_t>(uint64_t value) {
    this->value.i = value;
}
template<>
inline void Parameters::Parameter::set_internal<uint32_t>(uint32_t value) {
    this->value.i = value;
}
template<>
inline void Parameters::Parameter::set_internal<double>(double value) {
    this->value.d = value;
}
template<>
inline void Parameters::Parameter::set_internal<float>(float value) {
    this->value.d = value;
}
template<>
inline void Parameters::Parameter::set_internal<bool>(bool value) {
    this->value.b = value;
}
template<>
inline void Parameters::Parameter::set_internal<std::string>(std::string value) {
    this->value.s = new std::string(std::move(value));
}
template<>
inline void Parameters::Parameter::set_internal<const char*>(const char* value) {
    this->value.s = new std::string(value);
}
template<>
inline void Parameters::Parameter::set_internal<std::array<double, 3>>(std::array<double, 3> value) {
    this->value.da = value;
}
template<>
inline void Parameters::Parameter::set_internal<std::array<int64_t, 3>>(std::array<int64_t, 3> value) {
    this->value.ia = value;
}
template<>
inline void Parameters::Parameter::set_internal<ScoreFunction::AnchorScore>(ScoreFunction::AnchorScore value) {
    this->value.i = (int64_t) value;
}
template<>
inline void Parameters::Parameter::set_internal<logging::LoggingLevel>(logging::LoggingLevel value) {
    this->value.i = (int64_t) value;
}
template<>
inline void Parameters::Parameter::set_internal<Anchorer::ChainAlgorithm>(Anchorer::ChainAlgorithm value) {
    this->value.i = (int64_t) value;
}
template<>
inline void Parameters::Parameter::set_internal<Partitioner::ConstraintMethod>(Partitioner::ConstraintMethod value) {
    this->value.i = (int64_t) value;
}

template<typename T>
T Parameters::Parameter::get() const {
    check_type<T>();
    return get_internal<T>();
}

// default
template<typename T>
T Parameters::Parameter::get_internal() const {
    throw std::runtime_error(std::string("Unknown parameter type ") + typeid(T).name() +  " provided when getting parameter " + name);
    return T();
}
// specializations
template<>
inline int64_t Parameters::Parameter::get_internal<int64_t>() const {
    return value.i;
}
template<>
inline ScoreFunction::AnchorScore Parameters::Parameter::get_internal<ScoreFunction::AnchorScore>() const {
    return (ScoreFunction::AnchorScore) value.i;
}
template<>
inline logging::LoggingLevel Parameters::Parameter::get_internal<logging::LoggingLevel>() const {
    return (logging::LoggingLevel) value.i;
}
template<>
inline Anchorer::ChainAlgorithm Parameters::Parameter::get_internal<Anchorer::ChainAlgorithm>() const {
    return (Anchorer::ChainAlgorithm) value.i;
}
template<>
inline Partitioner::ConstraintMethod Parameters::Parameter::get_internal<Partitioner::ConstraintMethod>() const {
    return (Partitioner::ConstraintMethod) value.i;
}
template<>
inline double Parameters::Parameter::get_internal<double>() const {
    return value.d;
}
template<>
inline std::string Parameters::Parameter::get_internal<std::string>() const {
    return *value.s;
}
template<>
inline bool Parameters::Parameter::get_internal<bool>() const {
    return value.b;
}
template<>
inline std::array<int64_t, 3> Parameters::Parameter::get_internal<std::array<int64_t, 3>>() const {
    return value.ia;
}
template<>
inline std::array<double, 3> Parameters::Parameter::get_internal<std::array<double, 3>>() const {
    return value.da;
}


template<typename T>
void Parameters::enforce_lim(const std::string& name, bool lim_is_min, bool equal_ok, T lim) const {
    auto value = parameter(name).get<T>();
    if (!((value == lim && equal_ok) || (value > lim && lim_is_min) || (value < lim && !lim_is_min))) {
        std::stringstream strm;
        strm << "Invalid value for parameter '" << name << "', got " << value << " but expected a value ";
        if (lim_is_min) {
            strm << '>';
        }
        else {
            strm << '<';
        }
        if (equal_ok) {
            strm << '=';
        }
        strm << ' ' << lim;
        throw std::runtime_error(strm.str());
    }
}
template<typename T>
void Parameters::enforce_lt(const std::string& name, T upper_bound) const {
    enforce_lim(name, false, false, upper_bound);
}
template<typename T>
void Parameters::enforce_leq(const std::string& name, T upper_bound) const {
    enforce_lim(name, false, true, upper_bound);
}
template<typename T>
void Parameters::enforce_gt(const std::string& name, T lower_bound) const {
    enforce_lim(name, true, false, lower_bound);
}
template<typename T>
void Parameters::enforce_geq(const std::string& name, T lower_bound) const {
    enforce_lim(name, true, true, lower_bound);
}
template<typename T>
void Parameters::enforce_range(const std::string& name, T lower_bound, T upper_bound) const {
    enforce_geq(name, lower_bound);
    enforce_leq(name, upper_bound);
}

}

#endif /* centrolign_parameters_hpp */
