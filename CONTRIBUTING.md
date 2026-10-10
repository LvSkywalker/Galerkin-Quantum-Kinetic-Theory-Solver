# Contributing

Build and run CTest before submitting changes. For numerical modifications, describe the affected model, units and approximation, and add a check with a physical or numerical justification. Do not replace frozen references merely to make a failing test pass.

Changes to collision kernels should document their differential-cross-section convention and show conservation, scaling and relevant limiting cases. Changes to reported results should include the input settings, production outputs and regenerated figure data. Basis-size and quadrature studies should remain separate.

Use the repository's `.clang-format` for C++ changes. Keep generated binaries and local experiments under ignored `build*/` and `runs/` paths. When reporting a bug, include the command, compiler, platform, build type and relevant JSON diagnostics.

New methods and external data require source citations and license attribution. Performance claims should identify the hardware, compiler flags, problem size and timing procedure.
