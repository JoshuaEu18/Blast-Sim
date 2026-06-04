"""
Entry point for the Blast Damage Simulator.

Usage:
    cd blast_sim
    pip install -r requirements.txt
    python main.py

Then open  http://127.0.0.1:8050  in your browser.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from viz.dashboard import create_app

if __name__ == '__main__':
    app = create_app()
    print("\n" + "="*60)
    print("  Blast Damage Simulator")
    print("  Open:  http://127.0.0.1:8050")
    print("="*60 + "\n")
    app.run(debug=True, host='127.0.0.1', port=8050)
