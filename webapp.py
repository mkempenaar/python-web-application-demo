from flask import Flask, render_template, request, jsonify, make_response
from logic import translate_dna, generate_random_dna, mutate_dna

app = Flask(__name__)

## Template Rendering ##

@app.route('/')
def hello_world():
    """ Renders the DNA tools web form
    :return: the rendered view
    """
    return render_template('dnatools.html', title='Translate DNA Sequence')


@app.route('/dna-translate', methods=['POST'])
def process_webform():
    """ Receives a DNA sequence and selected frames and performs translation.
    The resulting dictionary with 'frame:translated sequence' is passed to the
    template and rendered as a table.
    :return: the rendered view containing the result table
    """
    # Retrieve the entered DNA-sequence (<textarea>)
    sequence = request.form.get('seq')

    # Retrieve the requested frames (<checkbox>; both forward and reverse)
    # and cast each list element to integer values
    fframes = [int(frame) for frame in request.form.getlist('fframe')]
    rframes = [int(frame) for frame in request.form.getlist('rframe')]

    # Translate the sequence with the requested frames
    translated_forward_seqs = translate_dna(sequence, fframes)
    translated_reverse_seqs = translate_dna(sequence, rframes, reverse=True)

    # Combine dictionaries
    translated_seqs = {**translated_forward_seqs, **translated_reverse_seqs}

    # Render result page
    return render_template('translated.html',
                           title="Translated DNA Sequence",
                           translated=translated_seqs)


## Data Providers ##

@app.route('/generate-dna', method=['GET'])
def generate_dna():
    """ Generates random DNA sequence(s) with a given length 
    input: number of sequences and the desired length (default 5 and 50)
    :return: the sequences in a list as JSON """
    nseq = request.args.get('nseq')
    length = request.args.get('length')

    sequences = generate_random_dna(nseq, length)
    return make_response(jsonify(sequences), 200)

@app.route('/mutate-dna', method=['GET'])
def mutate_dna():
    """ Mutate DNA sequence(s) given a probability that each base is mutated
    input: probability in percentage
    :return: possibly mutated sequence(s) in a list as JSON """
    prob = request.args.get('prob')
    sequences = request.args.get('sequences')

    mutated = mutate_dna(sequences, prob)
    return make_response(jsonify(mutated), 200)

## TODO: add functionality
@app.route('/some-path', method=['GET'])
def do_something():
    pass


if __name__ == '__main__':
    app.run()
