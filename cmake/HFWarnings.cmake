# hf_warnings: the warning set, as an INTERFACE target that Heavy Feather's own targets link
# PRIVATE, so none of it travels to a host that links the library.
#
# -Wconversion is on from the start: the prototype paid for every double -> int and
# size_t -> int by hand anyway, and it is far cheaper to keep a tree clean than to clean one.

add_library(hf_warnings INTERFACE)

if(HF_WARNINGS)
  if(MSVC)
    target_compile_options(hf_warnings INTERFACE /W4 /permissive-)
    if(HF_WARNINGS_AS_ERRORS)
      target_compile_options(hf_warnings INTERFACE /WX)
    endif()
  else()
    target_compile_options(hf_warnings INTERFACE
      -Wall -Wextra -Wpedantic -Wshadow -Wconversion
      -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual)
    if(HF_WARNINGS_AS_ERRORS)
      target_compile_options(hf_warnings INTERFACE -Werror)
    endif()
  endif()
endif()
