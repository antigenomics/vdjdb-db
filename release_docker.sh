build --no-cache -t vdjdbdb . 2>&1 | tee database/docker_build.log
run -v `pwd`/database/ vdjdbdb 2>&1 | tee database/docker_run.log