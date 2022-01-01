# Installing pspy

In order to build and run the Python Parasolid wrapper pspy,
this directory must contain the Parasolid library files from
your Parasolid distribution.

These should be named

`pskernel_archive_win_x64.lib`, `pskernel_archive_linux_x64.lib`, or `pskernel_archive_intel_macos.lib`

depending on your system. They can also be archived in zip directories to assist in committing to a 
repository; in this case, replace `.lib` with `.zip` in the above names, and the build files,
either setup.py or CMakeLists.txt, will extract the appropriate version for your system.