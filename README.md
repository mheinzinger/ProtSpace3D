# ProtSpace3D
This repository holds example code on how to visualize a set of proteins in 3D space. Proteins are first represented as high-dimensional vector representations (embeddings) derived from protein language models (pLMs) which are later projected to 3D for further analysis.
When hovering over points, images of structures predicted via AlphaFold 2 (AF2) are shown. An interactive 3D visualization of the structure can also be shown in the right panel via clicking on one of the points.

Using the family of Three-finger toxins (3FTx), we exemplify how such analysis can help to derive hypothesis on the phylogenetics of this family.
An example of this is hosted [here](http://3ftx.pythonanywhere.com/).

<br/>
<p align="center">
    <img width="70%" src="https://raw.githubusercontent.com/mheinzinger/ProtSpace3D/main/example_output.png" alt="Example output">
</p>
<br/>


# Getting started

Clone this repository and get started as described in the Usage section below:

```sh
git clone https://github.com/mheinzinger/ProtSpace3D.git
```

Set up a local python virtual environment and install dependencies:

```sh
python3 -m venv venv
. venv/bin/activate
pip install --upgrade pip wheel
pip install -r requirements.txt
```


# Usage

Simply run the dash board via:
```sh
python app.py
```
Next, you can access the 3D visulization by accessing http://127.0.0.1:8050/ in your browser.

