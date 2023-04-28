#include <cstdio>
#include <cstdlib>
#include <cstdlib>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <random>

#include "centrolign/tree.hpp"

using namespace std;
using namespace centrolign;




int main(int argc, char* argv[]) {
     
    random_device rd;
    default_random_engine gen(rd());
    
    // TODO: only testing that these don't crash, but the results
    // have been verified by eye using the debug output
    
    {
        Tree tree("((A,B),(C,D));");
    }

    {
        Tree tree("((A),(B,C));");
    }

    vector<string> wiki_examples{
        "(,,(,));",
        "(A,B,(C,D));",
        "(A,B,(C,D)E)F;",
        "(:0.1,:0.2,(:0.3,:0.4):0.5);",
        "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;",
        "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",
        "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
        "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;"
    };
    
    for (auto& newick_str : wiki_examples) {
        Tree tree(newick_str);
    }
    
    cerr << "passed all tests!" << endl;
}
