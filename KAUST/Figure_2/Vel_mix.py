#!/usr/bin/env python3

from PIL import Image

images = [
    Image.open("plt/nfields_00000.png"),
    Image.open("plt/nfields_00078.png"),
    Image.open("plt/nfields_00090.png"),
    Image.open("plt/nfields_00100.png"),
    Image.open("vel/nvel_00000.png"),
    Image.open("vel/nvel_00078.png"),
    Image.open("vel/nvel_00090.png"),
    Image.open("vel/nvel_00100.png"),
    ]

for i in range(0, 1):
    img = images[i]
    images[i] = img.crop((
        2,              # left
        0,            # top
        img.width,      # right
        img.height      # bottom
    ))
for i in range(4, 5):
    img = images[i]
    images[i] = img.crop((
        2,              # left
        0,            # top
        img.width,      # right
        img.height      # bottom
    ))

for i in range(4, 8):
    img = images[i]
    images[i] = img.crop((
        0,              # left
        200,            # top
        img.width,      # right
        img.height      # bottom
    ))

for i in range(0, 8):
    img = images[i]
    images[i] = img.crop((
        11.4,              # left
        0,            # top
        img.width-11.4,      # right
        img.height      # bottom
    ))

# Split into two rows
row1 = images[:4]
row2 = images[4:]

# Width of each column
col_widths = [
    max(row1[0].width, row2[0].width),
    max(row1[1].width, row2[1].width),
    max(row1[2].width, row2[2].width),
    max(row1[3].width, row2[3].width)
]

# Height of each row
row1_height = max(img.height for img in row1)
row2_height = max(img.height for img in row2)

total_width = sum(col_widths)
total_height = row1_height + row2_height
print(total_height,total_width)
combined = Image.new(
    "RGB",
    (total_width, total_height),
    "white"
)

# First row
x = 0
for i, img in enumerate(row1):
    combined.paste(img, (x, 0))
    x += col_widths[i]

# Second row
x = 0
for i, img in enumerate(row2):
    combined.paste(img, (x, row1_height))
    x += col_widths[i]

combined.save("combined.png", dpi=(600, 600))
combined.save("combined.pdf", dpi=(600, 600))
