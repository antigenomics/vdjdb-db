docker build --no-cache -t vdjdbdb . 2>&1 | tee database/docker_build.log
docker run -v `pwd`/database:/root/output vdjdbdb 2>&1 | tee database/docker_run.log