build:
	cp ../DESCRIPTION . && \
		docker build -t uscbiostats/netdiffuser:latest -f Dockerfile .
push: build
	docker push uscbiostats/netdiffuser:latest 
