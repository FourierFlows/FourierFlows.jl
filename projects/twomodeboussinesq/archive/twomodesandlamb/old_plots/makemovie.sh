#!/bin/bash

framename="$1_%03d.png"
filename="$1.mp4"

ffmpeg \
    -r 2 \
    -f image2 \
    -s 1920x1080 \
    -i $framename \
    -vcodec libx264 \
    -crf 25 \
    -pix_fmt yuv420p \
    $filename
