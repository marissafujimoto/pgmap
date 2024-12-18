import argparse
import sys


def _parse_args(args: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="pgmap", description="A tool to count paired guides from CRISPR double knockout screens.")
    parser.add_argument('-c', '--count', required=False,
                        type=_check_nonnegative)
    return parser.parse_args(args)


def _check_nonnegative(value: str) -> int:
    int_value = int(value)

    if int_value < 0:
        raise argparse.ArgumentTypeError(
            f"Count must be nonnegative but was {value}")

    return int_value


if __name__ == "__main__":
    _parse_args(sys.argv[1:])
