~/bin/google-compiler --js=lib/StarJs{.js,/Math.js,/Vector.js,/Time.js,/Coord.js,/Kepler.js,/Solar.js} --js=lib/export.js --js_output_file=StarJs.min.js --compilation_level ADVANCED_OPTIMIZATIONS
cat lib/StarJs{.js,/Math.js,/Vector.js,/Time.js,/Coord.js,/Kepler.js,/Solar.js} > StarJs.concat.js
