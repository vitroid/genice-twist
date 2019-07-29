.DELETE_ON_ERROR:
all: README.md

%: temp_% replacer.py genice_twist/formats/twist.py genice_twist/__init__.py
	python replacer.py < $< > $@


test: CRN1.btwc.test CRN1.CM.twist.png CRN1.DB.twist.png CRN1.SB.twist.png CRN1.twist.png \
      CRN1.twist.yap # CRN1.CM.twist.svg CRN1.DB.twist.svg CRN1.SB.twist.svg CRN1.twist.svg

%.test:
	make $*
	diff $* ref/$*
%.btwc: genice_twist/formats/twist.py
	genice $* -f twist > $@
%.CM.twist.png: genice_twist/formats/twist.py
	genice $* -f twist[png:CM] > $@
%.CM.twist.svg: genice_twist/formats/twist.py
	genice $* -f twist[svg:CM] > $@
%.DB.twist.png: genice_twist/formats/twist.py
	genice $* -f twist[png:DB] > $@
%.DB.twist.svg: genice_twist/formats/twist.py
	genice $* -f twist[svg:DB] > $@
%.SB.twist.png: genice_twist/formats/twist.py
	genice $* -f twist[png:SB] > $@
%.SB.twist.svg: genice_twist/formats/twist.py
	genice $* -f twist[svg:SB] > $@
%.twist.yap: genice_twist/formats/twist.py
	genice $* -f twist[yaplot] > $@
%.twist.png: genice_twist/formats/twist.py
	genice $* -f twist[png] > $@
%.twist.svg: genice_twist/formats/twist.py
	genice $* -f twist[svg] > $@

prepare: # might require root privilege.
	pip install genice genice-svg twist-op


test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install pairlist
	pip install --index-url https://test.pypi.org/simple/ genice-twist



install:
	./setup.py install
uninstall:
	-pip uninstall -y genice-twist
build: README.md $(wildcard genice_twist/formats/*.py)
	./setup.py sdist bdist_wheel


deploy: build
	twine upload dist/*
check:
	./setup.py check
clean:
	-rm $(ALL) *~ */*~ *.btwc *.png *.svg *.yap
	-rm -rf build dist *.egg-info
	-find . -name __pycache__ | xargs rm -rf
