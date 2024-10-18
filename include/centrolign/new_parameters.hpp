#ifndef centrolign_new_parameters_hpp
#define centrolign_new_parameters_hpp

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
struct NewParameters {
public:
    NewParameters();
//    NewParameters(std::istream& in);
    ~NewParameters() = default;
    
    // generate a config YAML that would recreate this set of parameters
    std::string generate_config() const;
    
//    // check logical validity of parameters, throw an error if anything is off
//    void validate() const;
//
//    // pass all of the parameters to configurable modules
//    void apply(Core& core) const;
    
    
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
    
    struct Parameter {
    public:
        
        template<typename T>
        Parameter(type_t type, const std::string& name, const std::string& help, T value);
        Parameter() = default;
        Parameter(Parameter&& other);
        Parameter(const Parameter& other);
        ~Parameter();
        
        Parameter& operator=(Parameter&& other) = delete;
        Parameter& operator=(const Parameter& other);
        
        template<typename T>
        void set(T value);
        
        type_t get_type() const;
        const std::string& get_name() const;
        const std::string& get_help() const;
        std::string value_str() const;
        
    private:
        
        Parameter(const std::string& name, const std::string& help) : name(name), help(help) { }
        
        template<typename T>
        void check_type() const;
        
        template<typename T>
        void set_internal(T val);
        
        const type_t type = (type_t) -1;
        const std::string name;
        const std::string help;
        value_t value = { 0 };
        
    };
        
    // the parameters, organized by 
    std::map<submodule_t, std::pair<std::string, std::vector<Parameter>>> params;
    std::unordered_map<std::string, std::pair<submodule_t, size_t>> param_position;
    
    void initialize_submodule(submodule_t module, const std::string& description);
    
    template<typename T>
    void add_parameter(submodule_t module, const std::string& name, type_t type, T default_value, const std::string& help);
    
    
    
};






/*
 * Template implementations
 */

template<typename T>
void NewParameters::add_parameter(submodule_t module, const std::string& name, type_t type, T default_value, const std::string& help) {
    
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
NewParameters::Parameter::Parameter(type_t type, const std::string& name, const std::string& help, T value) : type(type), name(name), help(help) {
    set<T>(value);
}

template<typename T>
void NewParameters::Parameter::check_type() const {
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
        (type == Enum && !std::is_enum<T>::value) ||
        (type == String && typeid(T) != typeid(std::string)) ||
        (type == DoubleArray3 && std::is_array<T>::value) ||
        (type == IntegerArray3 && std::is_array<T>::value)) {
        throw std::runtime_error(std::string("Value of type ") + typeid(value).name() + " for parameter " + name + " does not match expected type " + type_names.at(type));
    }
}

template<typename T>
void NewParameters::Parameter::set(T value) {
    check_type<T>();
    set_internal<T>(value);
}

// default
template<typename T>
void NewParameters::Parameter::set_internal(T value) {
    throw std::runtime_error(std::string("Parameter setting not implemented for type ") + typeid(value).name());
}
// specializations
template<>
void NewParameters::Parameter::set_internal<int64_t>(int64_t value) {
    this->value.i = value;
}
template<>
void NewParameters::Parameter::set_internal<int32_t>(int32_t value) {
    this->value.i = value;
}
template<>
void NewParameters::Parameter::set_internal<uint64_t>(uint64_t value) {
    this->value.i = value;
}
template<>
void NewParameters::Parameter::set_internal<uint32_t>(uint32_t value) {
    this->value.i = value;
}
template<>
void NewParameters::Parameter::set_internal<double>(double value) {
    this->value.d = value;
}
template<>
void NewParameters::Parameter::set_internal<float>(float value) {
    this->value.d = value;
}
template<>
void NewParameters::Parameter::set_internal<bool>(bool value) {
    this->value.b = value;
}
template<>
void NewParameters::Parameter::set_internal<std::string>(std::string value) {
    this->value.s = new std::string(std::move(value));
}
template<>
void NewParameters::Parameter::set_internal<std::array<double, 3>>(std::array<double, 3> value) {
    this->value.da = value;
}
template<>
void NewParameters::Parameter::set_internal<std::array<int64_t, 3>>(std::array<int64_t, 3> value) {
    this->value.ia = value;
}
template<>
void NewParameters::Parameter::set_internal<ScoreFunction::AnchorScore>(ScoreFunction::AnchorScore value) {
    this->value.i = (int64_t) value;
}
template<>
void NewParameters::Parameter::set_internal<logging::LoggingLevel>(logging::LoggingLevel value) {
    this->value.i = (int64_t) value;
}
template<>
void NewParameters::Parameter::set_internal<Anchorer::ChainAlgorithm>(Anchorer::ChainAlgorithm value) {
    this->value.i = (int64_t) value;
}
template<>
void NewParameters::Parameter::set_internal<Partitioner::ConstraintMethod>(Partitioner::ConstraintMethod value) {
    this->value.i = (int64_t) value;
}

}

#endif /* centrolign_new_parameters_hpp */
