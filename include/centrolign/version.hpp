#ifndef centrolign_version_hpp
#define centrolign_version_hpp

#include <string>

namespace centrolign {


// adapted from https://stackoverflow.com/questions/1435953/how-can-i-pass-git-sha1-to-compiler-as-definition-using-cmake

/*
 * Struct namespace to track the version of this software
 */
struct Version
{
    // major, minor, and patch version number for this release
    static const int MAJOR;
    static const int MINOR;
    static const int PATCH;
    // the git commit and info of the built version
    static const std::string GIT_HASH;
    static const std::string GIT_DATE;
    static const std::string GIT_COMMIT_SUBJECT;
};

}

#endif /* centrolign_version_hpp */
