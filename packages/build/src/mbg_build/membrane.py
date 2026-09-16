"""Membrane-building functions — port your COBY-based scripts here."""

import COBY

def build_membrane(*args, **kwargs):
    raise NotImplementedError("port your existing build script here")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    args = p.parse_args()
    build_membrane()
