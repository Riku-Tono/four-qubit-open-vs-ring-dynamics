"""Finite part of the noncancellation proof at lambda=1/2.

The report covers d >= 400 analytically and treats triangles, adjacent
pairs, and the four-site antipodes separately. Here all calculations are
integers after removing common positive channel-dependent factors.
"""


CUTOFF = 400


def binomial_row(length):
    row = [1]
    for k in range(1, length + 1):
        row.append(row[-1] * (length - k + 1) // k)
    return row


def signed_kernel(length, row, changes, include_yy=False):
    """Sum the distance-squared kernel on one parity of a binomial row."""
    total = 2**(length - 1)
    first = 0
    second = length * 2**(length - 3)
    fourth = (3 * length**2 - 2 * length) * 2**(length - 5) if include_yy else 0
    for position, decrement in changes:
        removed = decrement * row[position]
        offset = position - length // 2
        offset_squared = offset**2
        total -= removed
        first -= removed * offset
        second -= removed * offset_squared
        if include_yy:
            fourth -= removed * offset_squared**2
    zz = total * second - first**2
    if include_yy:
        yy = (length - 1)**2 * zz - (total * fourth - second**2)
        return yy, zz
    return zz


def path_scores(size, start, distance, xx_row, yz_row):
    """Signed scores of the clockwise arc, including its exterior bonds."""
    xx_changes, yz_changes = [], []
    xx_factor = yz_factor = 4
    # Edge k joins k to k+1 modulo N. Only the two edges at site 0
    # change the signed reciprocal weights: strong +1 -> -1, weak +1 -> -2.
    for edge, decrement in ((0, 2), (size - 1, 3)):
        offset = (edge - start) % size
        if offset < distance:
            xx_changes.append((2 * offset + 1, decrement))
            yz_changes.append((2 * offset + 2, decrement))
            if decrement == 3:
                xx_factor = yz_factor = 1
        elif offset in (distance, size - 1):
            position = 0 if offset == size - 1 else 2 * distance + 2
            yz_changes.append((position, decrement))
            if decrement == 3:
                yz_factor = 2

    if 0 in ((start - 1) % size, (start + distance + 1) % size):
        yz_factor = -yz_factor
    xx = signed_kernel(2 * distance, xx_row, xx_changes)
    yy, zz = signed_kernel(2 * distance + 2, yz_row, yz_changes, include_yy=True)
    return xx_factor * xx, yz_factor * yy, yz_factor * zz


def main():
    unique, tied = 0, 0
    for distance in range(2, CUTOFF):
        xx_row = binomial_row(2 * distance)
        yz_row = binomial_row(2 * distance + 2)

        size = 2 * distance + 1
        for start in range(size):
            scores = path_scores(size, start, distance, xx_row, yz_row)
            assert all(scores), (size, start, distance, scores)
            unique += 1

        if distance >= 3:
            size = 2 * distance
            for start in range(distance):
                first = path_scores(size, start, distance, xx_row, yz_row)
                second = path_scores(size, start + distance, distance, xx_row, yz_row)
                scores = tuple(a + b for a, b in zip(first, second))
                assert all(scores), (size, start, distance, scores)
                tied += 1

    print(f"No zero among {3 * unique:,} unique-path and {3 * tied:,} tied-path integer scores.")


if __name__ == "__main__":
    main()
