find . -type f ! -name '*.sh' -delete
rm -r -- ./*/
rm -r ../install/*
touch ../install/placeholder
