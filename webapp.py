#!/usr/bin/env python3

"""
Controller for a Python Flask and Bokeh demo project
"""
import os
import json
from bokeh.embed import json_item
from werkzeug.utils import secure_filename
from flask import Flask, render_template, request, jsonify, make_response
from graphics import create_gc_plot
from logic import translate_dna, generate_random_dna, mutate_dna
from logic import get_fasta_stats, calculate_molecular_weight

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = '/tmp'

## Template Rendering ##


@app.route('/')
def dnatools():
    """ Renders the DNA tools web form
    :return: the rendered view
    """
    return render_template('dnatools.html', title='Translate DNA Sequence')


## Data Providers ##

@app.route('/dna-translate-data', methods=['POST'])
def dna_translate_json():
    """ Receives a DNA sequence and selected frames and performs translation.
    Resulting translation is returned as a JSON string for further processing on the client.
    :return: JSON list containing a '(frame, protein sequence)' tuple for each selected frame """
    # Retrieve the entered DNA-sequence (<textarea>)
    sequence = request.form.get('seq')

    # Retrieve the requested frames (<checkbox>; both forward and reverse)
    # and cast each list element to integer values
    fframes = [int(frame) for frame in request.form.getlist('fframe')]
    rframes = [int(frame) for frame in request.form.getlist('rframe')]

    # Translate the sequence with the requested frames
    translated_forward_seqs = translate_dna(sequence, fframes)
    translated_reverse_seqs = translate_dna(sequence, rframes, reverse=True)

    # Combine frames
    translated_seqs = translated_forward_seqs + translated_reverse_seqs

    return make_response(jsonify(translated_seqs, 200))


@app.route('/generate-dna', methods=['GET'])
def generate_dna():
    """ Generates random DNA sequence(s) with a given length
    input: number of sequences and the desired length (default 5 and 60)
    If 'coding' is true, sequence is an open reading frame
    :return: the sequences in a list as JSON
    usage: http://127.0.0.1/generate-dna?nseq=5&length=100&coding=true
    """
    nseq = request.args.get('nseq', default=5)
    length = request.args.get('length', default=60)
    coding = request.args.get('coding', default=False)

    sequences = generate_random_dna(nseq, length, coding)
    return make_response(jsonify(sequences), 200)


@app.route('/mutate-dna', methods=['GET'])
def mutate_dna_sequence():
    """ Mutate DNA sequence(s) given a probability that each base is mutated
    input: single sequence and probability
    :return: possibly mutated sequence(s) in a list as JSON
    usage: http://127.0.0.1:5000/mutate-dna?prob=0.05&sequence=GTGAGAGAAAAGGCAGAGCTGGGCCAAGGCCCTGCCTCTCCGGGATGGTCTGTGGGGGAGCTGCAGCAGGGAGTG
    """
    prob = request.args.get('prob')
    sequence = request.args.get('sequence')

    mutated = mutate_dna(sequence, prob)
    return make_response(jsonify(mutated), 200)


@app.route('/fasta-statistics', methods=['POST'])
def fasta_statistics():
    """ Calculates basic sequence statistics for each sequence uploaded in a FASTA file:
    - GC-percentage
    - Molecular Weight
    - Sequence length
    - ...
    usage (using curl):
        curl -F 'file=@/full/path/to/data/sequences.fa' http://127.0.0.1:5000/fasta-statistics
    """
    filepath = save_uploaded_file(request, 'file')
    if not filepath:
        return make_response("No FASTA file given", 400)

    # Process file and get the statistics
    stats = get_fasta_stats(filepath)
    return make_response(jsonify(stats), 200)


@app.route('/molecular-weight', methods=['GET'])
def molecular_weight():
    """ Calculates the molecular weight of a single-strand DNA sequence """
    sequence = request.args.get('sequence')
    weight = calculate_molecular_weight(sequence)
    return make_response(weight, 200)


@app.route('/gc-percentage', methods=['POST'])
def gc_percentage():
    """ Creates a Bokeh plot showing GC-percentages for sequences included in
    a (multi-)FASTA file. """
    filepath = save_uploaded_file(request, 'file', ['.FASTA', '.FA'])
    if not filepath:
        return make_response("No FASTA file given", 400)

    plot = create_gc_plot(filepath)

    return json.dumps(json_item(plot))


## TODO: add functionality
## Hint: see https://www.bioinformatics.org/sms2/about.html
@app.route('/some-path', methods=['GET'])
def do_something():
    """ Rename the route and method name followed by adding logic """
    pass

## Utilities ##


def save_uploaded_file(request, form_field, supported_extensions):
    """ Saves an uploaded file, checking for supported file extensions """
    if request.method == 'POST':
        # check if the post request has the file part
        if form_field not in request.files:
            return False

    file = request.files[form_field]
    # Check if the uploaded file has the correct extension
    filename, file_extension = os.path.splitext(file.filename)
    if file_extension.upper() in supported_extensions:
        filename = secure_filename(file.filename)
        path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(path)
        return path
    
    return False


if __name__ == '__main__':
    app.secret_key = 'super secret key'
    app.run(debug=True)
