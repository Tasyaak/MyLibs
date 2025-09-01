#pragma once

#ifndef _MSVC_LANG
    #define _MSVC_LANG 0
#endif

// Позволяет отключить проверку при необходимости
#ifndef MYLIBS_DISABLE_STD_GUARD
    #if defined(_MSC_VER)
        #define MYLIBS_CPP_VER (_MSVC_LANG ? _MSVC_LANG : __cplusplus)
    #else
        #define MYLIBS_CPP_VER __cplusplus
#endif

    // Требуется C++17+
    static_assert(MYLIBS_CPP_VER >= 201703L,
        "MyLibs requires C++17 or later. Enable /std:c++17 (VS: C/C++ -> Language -> C++ Language Standard)");
#endif