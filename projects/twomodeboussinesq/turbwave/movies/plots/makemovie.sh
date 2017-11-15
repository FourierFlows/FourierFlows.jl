#!/bin/bash

framename="$1_%06d.png"
filename="$1.mp4"
framerate=8

ffmpeg \
    -r $framerate \
    -f image2 \
    -s 1920x1080 \
    -i $framename \
    -vcodec libx264 \
    -crf 25 \
    -pix_fmt yuv420p \
    $filename
