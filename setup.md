# The objective of this file is to present a fast way to setup a C++ financial instrument pricing environment on CLion JetBrains IDE.

## Setting up CLion : requirements ..

First starting by installing CLion JetBrains IDE from the following link : https://www.jetbrains.com/clion/download/#section=mac .
Make sure you install the right package whether you have intel-based or Apple Sillicon chip.
CLion includes an evaluation license key for a free 30-day trial.

The following link provides a quick tutorial on confirguring CLion on MacOS
https://www.jetbrains.com/help/clion/quick-tutorial-on-configuring-clion-on-macos.html#debugger .
Basically you will need Homebrew installed on your MacBook, which then you will use to install :
  - gcc compiler : "brew install gcc" should work, try "gcc --version" or "clang --version" to see which version of C++ you already have.
  - Make 
  -  ... check the website for more details.

The following link provides valuable tips on installing XCode, Homebrew and having a better understanding of the Mac Terminal.
https://www.moncefbelyamani.com/how-to-install-xcode-homebrew-git-rvm-ruby-on-mac/

See also : https://kig.re/2018/09/20/c++-newbie-tour-how-to-get-started-with-c++-on-mac-osx.html on how to get started with C++ for MacOSX.

### CMake and structuring projects 

In order to create good projects on C++, it is important to understand how projects are structured.
For this https://github.com/kigster/cmake-project-template provides a nice template to start with CMake.


## Setting QuantLib on CLion

This is also something I found for working with QuantLib on CLion : https://blog.jetbrains.com/clion/2015/12/quantlib-clion/ .
https://github.com/yhicc/FinancialLibrary .
Nice results in this google search : https://www.google.com/search?q=adding+quantlib+to+clion&rlz=1C5CHFA_enFR947FR947&oq=adding+quantlib+to+clion&aqs=chrome..69i57.3371j1j7&sourceid=chrome&ie=UTF-8 . 
But this is for later

https://github.com/yhicc/FinancialLibrary
