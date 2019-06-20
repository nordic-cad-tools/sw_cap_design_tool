from setuptools import setup

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name="sw_cap_design_tool",
    version="0.0.1",
    packages=["pyseeman"],
    url="",
    license="",
    author="thosou",
    author_email="",
    description="",
    install_requires=requirements,
)

