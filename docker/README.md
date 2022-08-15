To run the check in the main package folder:

```bash
docker run -v$(pwd):/fmcmc-pkg/ -w/fmcmc-pkg --rm -i uscbiostats/fmcmc-dev:latest make check
```

To update the version just type

```bash
make push
```
