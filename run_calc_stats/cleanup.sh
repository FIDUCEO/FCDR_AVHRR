
# find AVHRR*/ -name '*.log' | cpio -pdm logs/

find AVHRR*/ -name '*.log' -delete
find AVHRR*/ -name '*.sh' -delete
find AVHRR*/ -name '*.exe' -delete
find AVHRR*/ -name '*.py' -delete

find AVHRR*/ -name '*.data' -delete

